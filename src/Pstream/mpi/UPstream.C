/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "PstreamGlobals.H"
#include "profilingPstream.H"
#include "SubList.H"
#include "allReduce.H"
#include "int.H"
#include "collatedFileOperation.H"

#include <mpi.h>
// MediPack
#include <medi/medi.hpp>
#include <codi.hpp>
#include <codi/externals/codiMpiTypes.hpp>
using namespace medi;

#include <cstring>
#include <cstdlib>
#include <csignal>

#if defined(WM_SP)
    #define AMPI_SCALAR AMPI_FLOAT
    #define MPI_SOLVESCALAR MPI_FLOAT
#elif defined(WM_SPDP)
    #define AMPI_SCALAR AMPI_FLOAT
    #define MPI_SOLVESCALAR MPI_DOUBLE
#elif defined(WM_DP)
    #define AMPI_SCALAR AMPI_DOUBLE
    #define MPI_SOLVESCALAR MPI_DOUBLE
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// The min value and default for MPI buffers length
constexpr int minBufLen = 20000000;

// Track if we have attached MPI buffers
static bool ourBuffers = false;

// Track if we initialized MPI
static bool ourMpi = false;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

static void attachOurBuffers()
{
    if (ourBuffers)
    {
        return;  // Already attached
    }
    ourBuffers = true;

    // Use UPstream::mpiBufferSize (optimisationSwitch),
    // but allow override with MPI_BUFFER_SIZE env variable (int value)

#ifndef SGIMPI
    int len = 0;

    const std::string str(Foam::getEnv("MPI_BUFFER_SIZE"));
    if (str.empty() || !Foam::read(str, len) || len <= 0)
    {
        len = Foam::UPstream::mpiBufferSize;
    }

    if (len < minBufLen)
    {
        len = minBufLen;
    }

    if (Foam::UPstream::debug)
    {
        Foam::Pout<< "UPstream::init : buffer-size " << len << '\n';
    }

    char* buf = new char[len];

    if (MPI_SUCCESS != MPI_Buffer_attach(buf, len))
    {
        delete[] buf;
        Foam::Pout<< "UPstream::init : could not attach buffer\n";
    }
#endif
}


static void detachOurBuffers()
{
    if (!ourBuffers)
    {
        return;  // Nothing to detach
    }
    ourBuffers = false;

    // Some MPI notes suggest that the return code is MPI_SUCCESS when
    // no buffer is attached.
    // Be extra careful and require a non-zero size as well.

#ifndef SGIMPI
    int len = 0;
    char* buf = nullptr;

    if (MPI_SUCCESS == MPI_Buffer_detach(&buf, &len) && len)
    {
        delete[] buf;
    }
#endif
}

// whether to use Python, if yes, we do not call MPI_Finalize and let the 
// mpi4py finialize the MPI 
int isPython = 0;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// NOTE:
// valid parallel options vary between implementations, but flag common ones.
// if they are not removed by MPI_Init(), the subsequent argument processing
// will notice that they are wrong
void Foam::UPstream::addValidParOptions(HashTable<string>& validParOptions)
{
    validParOptions.insert("np", "");
    validParOptions.insert("p4pg", "PI file");
    validParOptions.insert("p4wd", "directory");
    validParOptions.insert("p4amslave", "");
    validParOptions.insert("p4yourname", "hostname");
    validParOptions.insert("machinefile", "machine file");
}


bool Foam::UPstream::initNull()
{
    int flag = 0;

    AMPI_Finalized(&flag);
    if (flag)
    {
        // Already finalized - this is an error
        FatalErrorInFunction
            << "MPI was already finalized - cannot perform MPI_Init\n"
            << Foam::abort(FatalError);

        return false;
    }

    AMPI_Initialized(&flag);
    if (flag)
    {
        // Already initialized - nothing to do
        // Initialize mpiTypes for AMPI_Datatype 
        PstreamGlobals::mpiTypes_ = new MpiTypes();
        return true;
    }
    else
    {
        // Not already initialized

    AMPI_Init_thread
    (
        nullptr,    // argc
        nullptr,    // argv
        AMPI_THREAD_SINGLE,
        &flag       // provided_thread_support
    );

    ourMpi = true;

    //AMPI_Init(nullptr, nullptr);

    // Initialize mpiTypes for AMPI_Datatype 
    PstreamGlobals::mpiTypes_ = new MpiTypes();

    return true;
}


bool Foam::UPstream::init(int& argc, char**& argv, const bool needsThread)
{
    int numprocs = 0, myRank = 0;
    int provided_thread_support = 0;
    // We need to check if the argv contains the -python option
    // if yes, we set isPython=1 and do not call MPI_Finalize, instead
    // we let mpi4py finalize the MPI
    // NOTE: this function is not called for serial runs, so we need some
    // special treatment in the ::exit function

    for(label i=0;i<argc;i++)
    {
        if(word(argv[i])=="-python") isPython=1;
    }

    int flag = 0;

    AMPI_Finalized(&flag);
    if (flag)
    {
        // Already finalized - this is an error
        FatalErrorInFunction
            << "MPI was already finalized - cannot perform MPI_Init" << endl
            << Foam::abort(FatalError);

        return false;
    }

    int provided_thread_support;

    AMPI_Initialized(&flag);
    if (flag)
    {
        // Already initialized - issue warning and skip the rest
        //WarningInFunction
        //    << "MPI was already initialized - cannot perform MPI_Init" << nl
        //    << "This could indicate an application programming error!" << endl;

        //return true;

        // NOTE: If MPI is initialized, call the AMPI_Init_common function
        // to initialize MeDiPack, check the AMPI_Init function defined
        // in MeDiPack/include/medi/ampi/wrappers.hpp
        AMPI_Init_common();
        // NOTE: Get the level of thread support provided. This is the same value that was 
        // returned in the provided argument in AMPI_Init_thread. 
        // provided_thread_support will be used later in setParRun
        AMPI_Query_thread(&provided_thread_support);
        // Initialize mpiTypes for AMPI_Datatype 
        PstreamGlobals::mpiTypes_ = new MpiTypes();
    }
    else
    {
        // If not initialized, do it here
        // AMPI_Init(&argc, &argv);
        AMPI_Init_thread
        (
            &argc,
            &argv,
            (
                needsThread
              ? AMPI_THREAD_MULTIPLE
              : AMPI_THREAD_SINGLE
            ),
            &provided_thread_support
        );
        // Initialize mpiTypes for AMPI_Datatype 
        PstreamGlobals::mpiTypes_ = new MpiTypes();
    }

    AMPI_Comm_size(AMPI_COMM_WORLD, &numprocs);
    AMPI_Comm_rank(AMPI_COMM_WORLD, &myRank);

    if (debug)
    {
        Pout<< "UPstream::init :"
            << " thread-support : wanted:" << needsThread
            << " obtained:"
            <<  (
                    provided_thread_support == MPI_THREAD_MULTIPLE
                  ? "MPI_THREAD_MULTIPLE"
                  : "MPI_THREAD_SINGLE"
                )
            << " procs:" << numprocs
            << " rank:" << myRank
            << " world:" << world << endl;
    }

    // Initialise parallel structure
    setParRun(numprocs, provided_thread_support == AMPI_THREAD_MULTIPLE);

    attachOurBuffers();

    return true;
}


void Foam::UPstream::shutdown(int errNo)
{
    // NOTE: return here for Python, we do not quit MPI in the OpenFOAM layer
    // NOTE: for serial runs, the ::init function will not be called so 
    // isPython will be always 0, we need some special treatment below
    if(isPython) return;

    if (debug)
    {
        Pout<< "UPstream::shutdown\n";
    }

    int flag = 0;

    AMPI_Initialized(&flag);
    if (!flag)
    {
        // No MPI initialized - we are done
        return;
    }
    else
    {
       // NOTE: here is the special treatment, if Python is used and if the
       // run is serial, the MPI_INIT WILL BE called by mpi4py, so we need
       // to see if the nProcs is 1, if yes, then return without finalizing
       // the MPI, again, we let mpi4py to finalize it
       int numprocs;
       AMPI_Comm_size(AMPI_COMM_WORLD, &numprocs);
       if(numprocs == 1) return; 
    }

    AMPI_Finalized(&flag);
    if (flag)
    {
        // Already finalized elsewhere?
        if (ourMpi)
        {
            WarningInFunction
                << "MPI was already finalized (by a connected program?)\n";
        }
        else if (debug)
        {
            Pout<< "UPstream::shutdown : was already finalized\n";
        }
    }
    else
    {
        detachOurBuffers();
    }


        flag = AMPI_Buffer_detach(&buf, &bufSize);

        if (AMPI_SUCCESS == flag && bufSize)
        {
            if (!PstreamGlobals::freedRequests_.found(requestID))
            {
                ++nOutstanding;
            }
        }

        PstreamGlobals::outstandingRequests_.clear();

        if (nOutstanding)
        {
            WarningInFunction
                << "There were still " << nOutstanding
                << " outstanding MPI_Requests." << nl
                << "Which means your code exited before doing a "
                << " UPstream::waitRequests()." << nl
                << "This should not happen for a normal code exit."
                << nl;
        }
    }

    // Clean mpi communicators
    forAll(myProcNo_, communicator)
    {
        if (myProcNo_[communicator] != -1)
        {
            freePstreamCommunicator(communicator);
        }
    }

    if (errnum == 0)
    {
        AMPI_Finalize();
        ::exit(errnum);
    }
    else
    {
        AMPI_Abort(AMPI_COMM_WORLD, errnum);
    }
}


void Foam::UPstream::exit(int errNo)
{
    UPstream::shutdown(errNo);
    std::exit(errNo);
}


void Foam::UPstream::abort()
{
    AMPI_Abort(AMPI_COMM_WORLD, 1);
}


void Foam::reduce
(
    scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    if (UPstream::parRun())
    {
        allReduce(Value, 1, PstreamGlobals::mpiTypes_->MPI_TYPE, AMPI_SUM, bop, tag, communicator);
    }
}


void Foam::reduce
(
    scalar& Value,
    const minOp<scalar>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    if (UPstream::parRun())
    {
        allReduce(Value, 1, PstreamGlobals::mpiTypes_->MPI_TYPE, AMPI_MIN, bop, tag, communicator);
    }
}


void Foam::reduce
(
    vector2D& Value,
    const sumOp<vector2D>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    if (UPstream::parRun())
    {
        allReduce(Value, 2, PstreamGlobals::mpiTypes_->MPI_TYPE, AMPI_SUM, bop, tag, communicator);
    }
}


void Foam::sumReduce
(
    scalar& Value,
    label& Count,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** sumReduce:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    vector2D twoScalars(Value, scalar(Count));
    reduce(twoScalars, sumOp<vector2D>(), tag, communicator);

    Value = twoScalars.x();
    Count = twoScalars.y().getValue();
}


void Foam::reduce
(
    scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label communicator,
    label& requestID
)
{
#ifdef MPIX_COMM_TYPE_SHARED
    // Assume mpich2 with non-blocking collectives extensions. Once mpi3
    // is available this will change.
    AMPI_Request request;
    scalar v = Value;
    AMPIX_Ireduce
    (
        &v,
        &Value,
        1,
        AMPI_SCALAR,
        AMPI_SUM,
        0,              //root
        PstreamGlobals::MPICommunicators_[communicator],
        &request
    );
}


#if defined(WM_SPDP)
void Foam::reduce
(
    solveScalar& Value,
    const sumOp<solveScalar>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    allReduce(Value, 1, MPI_SOLVESCALAR, MPI_SUM, bop, tag, communicator);
}


void Foam::reduce
(
    solveScalar& Value,
    const minOp<solveScalar>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    allReduce(Value, 1, MPI_SOLVESCALAR, MPI_MIN, bop, tag, communicator);
}


void Foam::reduce
(
    Vector2D<solveScalar>& Value,
    const sumOp<Vector2D<solveScalar>>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    allReduce(Value, 2, MPI_SOLVESCALAR, MPI_SUM, bop, tag, communicator);
}


void Foam::sumReduce
(
    solveScalar& Value,
    label& Count,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    Vector2D<solveScalar> twoScalars(Value, solveScalar(Count));
    reduce(twoScalars, sumOp<Vector2D<solveScalar>>(), tag, communicator);

    Value = twoScalars.x();
    Count = twoScalars.y();
}


void Foam::reduce
(
    solveScalar& Value,
    const sumOp<solveScalar>& bop,
    const int tag,
    const label communicator,
    label& requestID
)
{
    iallReduce<solveScalar>
    (
        &Value,
        1,
        MPI_SOLVESCALAR,
        MPI_SUM,
        communicator,
        requestID
    );
}


void Foam::reduce
(
    solveScalar values[],
    const int size,
    const sumOp<solveScalar>& bop,
    const int tag,
    const label communicator,
    label& requestID
)
{
    iallReduce<solveScalar>
    (
        values,
        size,
        MPI_SOLVESCALAR,
        MPI_SUM,
        communicator,
        requestID
    );
}
#endif


void Foam::UPstream::allToAll
(
    const labelUList& sendData,
    labelUList& recvData,
    const label communicator
)
{
    const label np = nProcs(communicator);

    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** allToAll :"
            << " np:" << np
            << " sendData:" << sendData.size()
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if (sendData.size() != np || recvData.size() != np)
    {
        FatalErrorInFunction
            << "Size of sendData " << sendData.size()
            << " or size of recvData " << recvData.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        recvData.deepCopy(sendData);
    }
    else
    {
        profilingPstream::beginTiming();

        if
        (
            // CoDiPack4OpenFOAM TODO Alltoall function is not AMPI yet
            // This shouldn't be an issue since the allToAll function is only used in
            // src/OpenFOAM/db/IOstreams/Pstreams/exchange.C to exchange sizes 
            MPI_Alltoall
            (
                // NOTE: const_cast is a temporary hack for
                // backward-compatibility with versions of OpenMPI < 1.7.4
                const_cast<label*>(sendData.cdata()),
                sizeof(label),
                MPI_BYTE,
                recvData.data(),
                sizeof(label),
                MPI_BYTE,
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Alltoall failed for " << sendData
                << " on communicator " << communicator
                << Foam::abort(FatalError);
        }

        profilingPstream::addAllToAllTime();
    }
}


void Foam::UPstream::allToAll
(
    const char* sendData,
    const UList<int>& sendSizes,
    const UList<int>& sendOffsets,
    char* recvData,
    const UList<int>& recvSizes,
    const UList<int>& recvOffsets,
    const word callerInfo,
    const std::type_info& typeInfo,
    const label communicator
)
{
    const label np = nProcs(communicator);

    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** MPI_Alltoallv :"
            << " sendSizes:" << sendSizes
            << " sendOffsets:" << sendOffsets
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if
    (
        sendSizes.size() != np
     || sendOffsets.size() != np
     || recvSizes.size() != np
     || recvOffsets.size() != np
    )
    {
        FatalErrorInFunction
            << "Size of sendSize " << sendSizes.size()
            << ", sendOffsets " << sendOffsets.size()
            << ", recvSizes " << recvSizes.size()
            << " or recvOffsets " << recvOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    bool typeActive = Foam::PstreamGlobals::isTypeActive(typeInfo)
                   && codi::RealReverse::getGlobalTape().isActive();

    if (debug)
    {
        Pout<< "UPstream::allToAll :"
            << " typeActive: " << typeActive << " typeid: " << typeInfo.name()
            << Foam::endl;
    }

    if (!UPstream::parRun())
    {
        if (recvSizes[0] != sendSizes[0])
        {
            FatalErrorInFunction
                << "Bytes to send " << sendSizes[0]
                << " does not equal bytes to receive " << recvSizes[0]
                << Foam::abort(FatalError);
        }
        std::memmove(recvData, &sendData[sendOffsets[0]], recvSizes[0]);
    }
    else
    {
        label Err = 0;
        if (typeActive)
        {
            Err = AMPI_Alltoallv
            (
                reinterpret_cast<scalar*>(const_cast<char*>(sendData)),
                const_cast<int*>(sendSizes.begin()),
                const_cast<int*>(sendOffsets.begin()),
                PstreamGlobals::mpiTypes_->MPI_TYPE,
                reinterpret_cast<scalar*>(recvData),
                const_cast<int*>(recvSizes.begin()),
                const_cast<int*>(recvOffsets.begin()),
                PstreamGlobals::mpiTypes_->MPI_TYPE,
                PstreamGlobals::MPICommunicators_[communicator]
            );
        }
        else
        {
            Err = AMPI_Alltoallv
            (
                reinterpret_cast<unsigned char*>(const_cast<char*>(sendData)),
                const_cast<int*>(sendSizes.begin()),
                const_cast<int*>(sendOffsets.begin()),
                AMPI_BYTE,
                reinterpret_cast<unsigned char*>(recvData),
                const_cast<int*>(recvSizes.begin()),
                const_cast<int*>(recvOffsets.begin()),
                AMPI_BYTE,
                PstreamGlobals::MPICommunicators_[communicator]
            );
        }
        if (Err)
        {
            FatalErrorInFunction
                << "MPI_Alltoallv failed for sendSizes " << sendSizes
                << " recvSizes " << recvSizes
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

        profilingPstream::addAllToAllTime();
    }
}


void Foam::UPstream::mpiGather
(
    const char* sendData,
    int sendSize,

    char* recvData,
    int recvSize,
    const label communicator
)
{
    const label np = nProcs(communicator);

    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** MPI_Gather :"
            << " np:" << np
            << " recvSize:" << recvSize
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if (!UPstream::parRun())
    {
        std::memmove(recvData, sendData, recvSize);
    }
    else
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Gather
            (
                const_cast<char*>(sendData),
                sendSize,
                MPI_BYTE,
                recvData,
                recvSize,
                MPI_BYTE,
                0,
                MPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Gather failed for sendSize " << sendSize
                << " recvSize " << recvSize
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

        profilingPstream::addGatherTime();
    }
}


void Foam::UPstream::mpiScatter
(
    const char* sendData,
    int sendSize,

    char* recvData,
    int recvSize,
    const label communicator
)
{
    const label np = nProcs(communicator);

    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** MPI_Scatter :"
            << " np:" << np
            << " recvSize:" << recvSize
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if (!UPstream::parRun())
    {
        std::memmove(recvData, sendData, recvSize);
    }
    else
    {
        profilingPstream::beginTiming();

        if
        (
            MPI_Scatter
            (
                const_cast<char*>(sendData),
                sendSize,
                MPI_BYTE,
                recvData,
                recvSize,
                MPI_BYTE,
                0,
                MPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Scatter failed for sendSize " << sendSize
                << " recvSize " << recvSize
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

        profilingPstream::addScatterTime();
    }
}


void Foam::UPstream::gather
(
    const char* sendData,
    int sendSize,

    char* recvData,
    const UList<int>& recvSizes,
    const UList<int>& recvOffsets,
    const word callerInfo,
    const std::type_info& typeInfo,
    const label communicator
)
{
    const label np = nProcs(communicator);

    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** MPI_Gatherv :"
            << " np:" << np
            << " recvSizes:" << recvSizes
            << " recvOffsets:" << recvOffsets
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if
    (
        UPstream::master(communicator)
     && (recvSizes.size() != np || recvOffsets.size() < np)
    )
    {
        // Note: allow recvOffsets to be e.g. 1 larger than np so we
        // can easily loop over the result

        FatalErrorInFunction
            << "Size of recvSizes " << recvSizes.size()
            << " or recvOffsets " << recvOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    bool typeActive = Foam::PstreamGlobals::isTypeActive(typeInfo)
                   && codi::RealReverse::getGlobalTape().isActive();

    if (debug)
    {
        Pout<< "UPstream::gather :"
            << " typeActive: " << typeActive << " typeid: " << typeInfo.name()
            << Foam::endl;
    }

    if (!UPstream::parRun())
    {
        // recvSizes[0] may be invalid - use sendSize instead
        std::memmove(recvData, sendData, sendSize);
    }
    else
    {
        label Err = 0;
        if (typeActive)
        {
            Err = AMPI_Gatherv
            (
                reinterpret_cast<scalar*>(const_cast<char*>(sendData)),
                sendSize,
                PstreamGlobals::mpiTypes_->MPI_TYPE,
                reinterpret_cast<scalar*>(recvData),
                const_cast<int*>(recvSizes.begin()),
                const_cast<int*>(recvOffsets.begin()),
                PstreamGlobals::mpiTypes_->MPI_TYPE,
                0,
                AMPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            );
        }
        else
        {
            Err = AMPI_Gatherv
            (
                reinterpret_cast<unsigned char*>(const_cast<char*>(sendData)),
                sendSize,
                AMPI_BYTE,
                reinterpret_cast<unsigned char*>(recvData),
                const_cast<int*>(recvSizes.begin()),
                const_cast<int*>(recvOffsets.begin()),
                AMPI_BYTE,
                0,
                AMPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            );
        }

        if (Err)
        {
            FatalErrorInFunction
                << "MPI_Gatherv failed for sendSize " << sendSize
                << " recvSizes " << recvSizes
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

        profilingPstream::addGatherTime();
    }
}


void Foam::UPstream::scatter
(
    const char* sendData,
    const UList<int>& sendSizes,
    const UList<int>& sendOffsets,

    char* recvData,
    int recvSize,
    const word callerInfo,
    const std::type_info& typeInfo,
    const label communicator
)
{
    const label np = nProcs(communicator);

    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** MPI_Scatterv :"
            << " np:" << np
            << " sendSizes:" << sendSizes
            << " sendOffsets:" << sendOffsets
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    if
    (
        UPstream::master(communicator)
     && (sendSizes.size() != np || sendOffsets.size() != np)
    )
    {
        FatalErrorInFunction
            << "Size of sendSizes " << sendSizes.size()
            << " or sendOffsets " << sendOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    bool typeActive = Foam::PstreamGlobals::isTypeActive(typeInfo)
                   && codi::RealReverse::getGlobalTape().isActive();

    if (debug)
    {
        Pout<< "UPstream::scatter :"
            << " typeActive: " << typeActive << " typeid: " << typeInfo.name()
            << Foam::endl;
    }

    if (!UPstream::parRun())
    {
        std::memmove(recvData, sendData, recvSize);
    }
    else
    {
        label Err = 0;
        if (typeActive)
        {
            Err = AMPI_Scatterv
            (
                reinterpret_cast<scalar*>(const_cast<char*>(sendData)),
                const_cast<int*>(sendSizes.begin()),
                const_cast<int*>(sendOffsets.begin()),
                PstreamGlobals::mpiTypes_->MPI_TYPE,
                reinterpret_cast<scalar*>(recvData),
                recvSize,
                PstreamGlobals::mpiTypes_->MPI_TYPE,
                0,
                AMPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            );
        }
        else
        {
            Err = AMPI_Scatterv
            (
                reinterpret_cast<unsigned char*>(const_cast<char*>(sendData)),
                const_cast<int*>(sendSizes.begin()),
                const_cast<int*>(sendOffsets.begin()),
                AMPI_BYTE,
                reinterpret_cast<unsigned char*>(recvData),
                recvSize,
                AMPI_BYTE,
                0,
                AMPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            );
        }

        if (Err)
        {
            FatalErrorInFunction
                << "MPI_Scatterv failed for sendSizes " << sendSizes
                << " sendOffsets " << sendOffsets
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

        profilingPstream::addScatterTime();
    }
}


void Foam::UPstream::allocatePstreamCommunicator
(
    const label parentIndex,
    const label index
)
{
    if (index == PstreamGlobals::MPIGroups_.size())
    {
        // Extend storage with dummy values
        AMPI_Group newGroup = AMPI_GROUP_NULL;
        PstreamGlobals::MPIGroups_.append(newGroup);
        AMPI_Comm newComm = AMPI_COMM_NULL;
        PstreamGlobals::MPICommunicators_.append(newComm);
    }
    else if (index > PstreamGlobals::MPIGroups_.size())
    {
        FatalErrorInFunction
            << "PstreamGlobals out of sync with UPstream data. Problem."
            << Foam::exit(FatalError);
    }


    if (parentIndex == -1)
    {
        // Allocate world communicator

        if (index != UPstream::worldComm)
        {
            FatalErrorInFunction
                << "world communicator should always be index "
                << UPstream::worldComm << Foam::exit(FatalError);
        }

        PstreamGlobals::MPICommunicators_[index] = AMPI_COMM_WORLD;
        AMPI_Comm_group(AMPI_COMM_WORLD, &PstreamGlobals::MPIGroups_[index]);
        AMPI_Comm_rank
        (
            PstreamGlobals::MPICommunicators_[index],
           &myProcNo_[index]
        );

        // Set the number of processes to the actual number
        int numProcs;
        AMPI_Comm_size(PstreamGlobals::MPICommunicators_[index], &numProcs);

        //procIDs_[index] = identity(numProcs);
        procIDs_[index].setSize(numProcs);
        forAll(procIDs_[index], i)
        {
            procIDs_[index][i] = i;
        }
    }
    else
    {
        // Create new group
        AMPI_Group_incl
        (
            PstreamGlobals::MPIGroups_[parentIndex],
            procIDs_[index].size(),
            procIDs_[index].cdata(),
           &PstreamGlobals::MPIGroups_[index]
        );

        // Create new communicator
        AMPI_Comm_create
        (
            PstreamGlobals::MPICommunicators_[parentIndex],
            PstreamGlobals::MPIGroups_[index],
            &PstreamGlobals::MPICommunicators_[index]
        );
        #else
        // Create new communicator for this group
        MPI_Comm_create_group
        (
            PstreamGlobals::MPICommunicators_[parentIndex],
            PstreamGlobals::MPIGroups_[index],
            Pstream::msgType(),
           &PstreamGlobals::MPICommunicators_[index]
        );
        #endif

        if (PstreamGlobals::MPICommunicators_[index] == AMPI_COMM_NULL)
        {
            myProcNo_[index] = -1;
        }
        else
        {
            if
            (
                AMPI_Comm_rank
                (
                    PstreamGlobals::MPICommunicators_[index],
                   &myProcNo_[index]
                )
            )
            {
                FatalErrorInFunction
                    << "Problem :"
                    << " when allocating communicator at " << index
                    << " from ranks " << procIDs_[index]
                    << " of parent " << parentIndex
                    << " cannot find my own rank"
                    << Foam::exit(FatalError);
            }
        }
    }
}


void Foam::UPstream::freePstreamCommunicator(const label communicator)
{
    if (communicator != 0)
    {
        if (PstreamGlobals::MPICommunicators_[communicator] != AMPI_COMM_NULL)
        {
            // Free communicator. Sets communicator to MPI_COMM_NULL
            AMPI_Comm_free(&PstreamGlobals::MPICommunicators_[communicator]);
        }
        if (PstreamGlobals::MPIGroups_[communicator] != AMPI_GROUP_NULL)
        {
            // Free greoup. Sets group to MPI_GROUP_NULL
            AMPI_Group_free(&PstreamGlobals::MPIGroups_[communicator]);
        }
    }
}


Foam::label Foam::UPstream::nRequests()
{
    return PstreamGlobals::outstandingRequests_.size();
}


void Foam::UPstream::resetRequests(const label i)
{
    if (i < PstreamGlobals::outstandingRequests_.size())
    {
        PstreamGlobals::outstandingRequests_.setSize(i);
    }
}


void Foam::UPstream::waitRequests(const label start)
{
    if (UPstream::debug)
    {
        Pout<< "UPstream::waitRequests : starting wait for "
            << PstreamGlobals::outstandingRequests_.size()-start
            << " outstanding requests starting at " << start << endl;
    }

    if (PstreamGlobals::outstandingRequests_.size())
    {
        SubList<AMPI_Request> waitRequests
        (
            PstreamGlobals::outstandingRequests_,
            PstreamGlobals::outstandingRequests_.size() - start,
            start
        );

        profilingPstream::beginTiming();

        if
        (
            AMPI_Waitall
            (
                waitRequests.size(),
                waitRequests.begin(),
                AMPI_STATUSES_IGNORE
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Waitall returned with error" << Foam::endl;
        }

        profilingPstream::addWaitTime();

        resetRequests(start);
    }

    if (debug)
    {
        Pout<< "UPstream::waitRequests : finished wait." << endl;
    }
}


void Foam::UPstream::waitRequest(const label i)
{
    if (debug)
    {
        Pout<< "UPstream::waitRequest : starting wait for request:" << i
            << endl;
    }

    if (i < 0 || i >= PstreamGlobals::outstandingRequests_.size())
    {
        FatalErrorInFunction
            << "There are " << PstreamGlobals::outstandingRequests_.size()
            << " outstanding send requests and you are asking for i=" << i
            << nl
            << "Maybe you are mixing blocking/non-blocking comms?"
            << Foam::abort(FatalError);
    }

    profilingPstream::beginTiming();

    if
    (
        AMPI_Wait
        (
           &PstreamGlobals::outstandingRequests_[i],
            AMPI_STATUS_IGNORE
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Wait returned with error" << Foam::endl;
    }

    profilingPstream::addWaitTime();
    // Push index onto free cache
    PstreamGlobals::freedRequests_.append(i);

    if (debug)
    {
        Pout<< "UPstream::waitRequest : finished wait for request:" << i
            << endl;
    }
}


bool Foam::UPstream::finishedRequest(const label i)
{
    if (debug)
    {
        Pout<< "UPstream::finishedRequest : checking request:" << i
            << endl;
    }

    if (i >= PstreamGlobals::outstandingRequests_.size())
    {
        FatalErrorInFunction
            << "There are " << PstreamGlobals::outstandingRequests_.size()
            << " outstanding send requests and you are asking for i=" << i
            << nl
            << "Maybe you are mixing blocking/non-blocking comms?"
            << Foam::abort(FatalError);
    }

    int flag;
    AMPI_Test
    (
       &PstreamGlobals::outstandingRequests_[i],
       &flag,
        AMPI_STATUS_IGNORE
    );

    if (debug)
    {
        Pout<< "UPstream::finishedRequest : finished request:" << i
            << endl;
    }

    return flag != 0;
}


int Foam::UPstream::allocateTag(const char* s)
{
    int tag;
    if (PstreamGlobals::freedTags_.size())
    {
        tag = PstreamGlobals::freedTags_.remove();
    }
    else
    {
        tag = PstreamGlobals::nTags_++;
    }

    if (debug)
    {
        //if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = 'X';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::allocateTag " << s
            << " : tag:" << tag
            << endl;
    }

    return tag;
}


int Foam::UPstream::allocateTag(const word& s)
{
    int tag;
    if (PstreamGlobals::freedTags_.size())
    {
        tag = PstreamGlobals::freedTags_.remove();
    }
    else
    {
        tag = PstreamGlobals::nTags_++;
    }

    if (debug)
    {
        //if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = 'X';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::allocateTag " << s
            << " : tag:" << tag
            << endl;
    }

    return tag;
}


void Foam::UPstream::freeTag(const char* s, const int tag)
{
    if (debug)
    {
        //if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = ' ';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::freeTag " << s << " tag:" << tag << endl;
    }
    PstreamGlobals::freedTags_.append(tag);
}


void Foam::UPstream::freeTag(const word& s, const int tag)
{
    if (debug)
    {
        //if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = ' ';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::freeTag " << s << " tag:" << tag << endl;
    }
    PstreamGlobals::freedTags_.append(tag);
}

#include <medi/medi.cpp>
// ************************************************************************* //
