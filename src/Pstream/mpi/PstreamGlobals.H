/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2022-2025 OpenCFD Ltd.
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

Namespace
    Foam::PstreamGlobals

Description
    Global functions and variables for working with parallel streams,
    but principally for MPI.

SourceFiles
    PstreamGlobals.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PstreamGlobals_H
#define Foam_PstreamGlobals_H

#include "DynamicList.H"
#include "FixedList.H"
#include "UPstream.H"  // For UPstream::Request
#include "openfoam_mpi.H"

// codi: MediPack for AD-aware MPI
#include <medi/medi.hpp>
#include <codi.hpp>
#include <codi/tools/mpi/codiMpiTypes.hpp>

using namespace medi;

#ifdef CODI_ADF
using MpiTypes = codi::CoDiMpiTypes<codi::RealForward>;
#endif
#ifdef CODI_ADR
using MpiTypes = codi::CoDiMpiTypes<codi::RealReverse>;
#endif

// codi: Global MediPack type pointer (defined in PstreamGlobals.C)
extern MpiTypes* mpiTypes;

// codi: Helper function to convert MPI_Op to AMPI_Op
inline AMPI_Op toAMPI_Op(MPI_Op op)
{
    if (op == MPI_SUM) return AMPI_SUM;
    if (op == MPI_MAX) return AMPI_MAX;
    if (op == MPI_MIN) return AMPI_MIN;

    FatalErrorInFunction
        << "Unsupported MPI_Op for AMPI operations. "
        << "Only MPI_SUM, MPI_MAX, and MPI_MIN are supported."
        << Foam::abort(Foam::FatalError);

    return AMPI_SUM;  // Never reached
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamGlobals
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Track if MPI_Comm_free is needed for communicator index in MPICommunicators_
extern DynamicList<bool> pendingMPIFree_;

// Current communicators, which may be allocated or predefined
// (eg, MPI_COMM_SELF, MPI_COMM_WORLD)
extern DynamicList<MPI_Comm> MPICommunicators_;

//- codi: Outstanding non-blocking operations (stores both MPI and AMPI requests).
extern DynamicList<AMPI_Request> outstandingRequests_;

// The fundamental count for each UPstream::dataTypes entry
typedef Foam::FixedList<int, 17> DataTypeCountLookupTable;

//- Fundamental count for each valid UPstream::dataTypes entry
//- Indexed by UPstream::dataTypes enum
extern DataTypeCountLookupTable dataTypesCount_;

// For UPstream::dataTypes lookup, includes space for last 'invalid' entry
typedef Foam::FixedList<MPI_Datatype, 17> DataTypeLookupTable;

//- MPI data types corresponding to fundamental and OpenFOAM types.
//- Indexed by UPstream::dataTypes enum
extern DataTypeLookupTable MPIdataTypes_;

// For UPstream::opCodes lookup, includes space for last 'invalid' entry
typedef Foam::FixedList<MPI_Op, 13> OpCodesLookupTable;

//- MPI operation types, indexed by UPstream::opCodes enum
extern OpCodesLookupTable MPIopCodes_;



// * * * * * * * * * * * * * * * Communicators * * * * * * * * * * * * * * * //

//- Initialize bookkeeping for MPI communicator index
void initCommunicator(const label index);

//- Fatal if communicator is outside the allocated range
inline void checkCommunicator(int comm, int rank)
{
    if (FOAM_UNLIKELY(comm < 0 || comm >= MPICommunicators_.size()))
    {
        FatalErrorInFunction
            << "rank:" << rank << " : illegal communicator "
            << comm << nl
            << "Communicator should be within range [0,"
            << MPICommunicators_.size()
            << ')' << Foam::abort(FatalError);
    }
}

//- True if warn communicator is active and not equal to given communicator
inline bool warnCommunicator(int comm) noexcept
{
    return (UPstream::warnComm >= 0 && comm != UPstream::warnComm);
}


// * * * * * * * * * * * * * * * * Data Types  * * * * * * * * * * * * * * * //

//- Create mapping into MPIdataTypes_ and define user data types
void initDataTypes();

//- Free any user data types
void deinitDataTypes();

//- Debugging only: check if data type mappings are non-null
bool checkDataTypes();

//- Debugging only: print data type names (all or just user-defined)
void printDataTypes(bool all = false);

//- Lookup of dataTypes enumeration as an MPI_Datatype
inline MPI_Datatype getDataType(UPstream::dataTypes id)
{
    return MPIdataTypes_[static_cast<int>(id)];
}

//- Fatal if data type is not valid
inline void checkDataType(UPstream::dataTypes id)
{
    if (id == UPstream::dataTypes::invalid)
    {
        FatalErrorInFunction
            << "Invalid data type"
            << Foam::abort(FatalError);
    }
}

//- Return MPI internal name for specified MPI_Datatype
std::string dataType_name(MPI_Datatype datatype);

//- Return MPI internal name for dataTypes enumeration
inline std::string dataType_name(UPstream::dataTypes id)
{
    return dataType_name(MPIdataTypes_[static_cast<int>(id)]);
}


// * * * * * * * * * * * * * * * * Op Codes  * * * * * * * * * * * * * * * * //

//- Create mapping into MPIopCodes_
void initOpCodes();

//- Free any user-defined op codes
void deinitOpCodes();

//- Debugging only: check if op code mappings are non-null
bool checkOpCodes();

//- Lookup of opCodes enumeration as an MPI_Op
inline MPI_Op getOpCode(UPstream::opCodes id)
{
    return MPIopCodes_[static_cast<int>(id)];
}

//- Fatal if opcode is not valid
inline void checkOpCode(UPstream::opCodes id)
{
    if (id == UPstream::opCodes::invalid)
    {
        FatalErrorInFunction
            << "Invalid operation code"
            << Foam::abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Requests  * * * * * * * * * * * * * * * * //

//- Reset UPstream::Request to MPI_REQUEST_NULL
//  Does not affect the stack of outstanding requests
inline void reset_request(UPstream::Request* req) noexcept
{
    if (req) *req = UPstream::Request(MPI_REQUEST_NULL);
}


// codi: the new push_request to take AMPI_Request as an input
inline void push_request
(
    AMPI_Request request,
    UPstream::Request* req = nullptr
)
{
    // here AMPI_Request is a struct AMPI_Request.request is an MPI_Request
    if (MPI_REQUEST_NULL != request.request)
    {
        // Always push full AMPI_Request to global list to preserve AD metadata
        label index = PstreamGlobals::outstandingRequests_.size();
        PstreamGlobals::outstandingRequests_.push_back(request);

        // NOTE: here we need special treatment if req is provided. 
        // we Store NEGATIVE index to distinguish from positive MPI_Request handles
        // Index 0 -> -1, index 1 -> -2, index 2 -> -3, etc.
        // so if the index is negative, it corresponds to AMPI_Request
        // if the index is positive, it corresponds to MPI_Request
        if (req)
        {
            *req = UPstream::Request(-(index + 1));
        }
    }
}

// codi: the original push_request that used MPI_Request, we will convert it to AMPI_Request
inline void push_request
(
    MPI_Request request,
    UPstream::Request* req = nullptr
)
{
    if (req)
    {
        // Request& API: Transcribe as UPstream::Request
        *req = UPstream::Request(request);
    }
    else if (MPI_REQUEST_NULL != request)
    {
        // Label-based API: Wrap MPI_Request in AMPI_Request and push to unified list
        // This ensures waitRequestAtIndex() can find it
        AMPI_Request ampiReq;
        ampiReq.request = request;
        ampiReq.handle = nullptr;
        ampiReq.func = nullptr;
        ampiReq.start = nullptr;
        ampiReq.end = nullptr;
        ampiReq.isActive = false;
        ampiReq.reverseData = nullptr;
        ampiReq.deleteDataFunc = nullptr;

        // Note: outstandingRequests_ saves only AMPI_Request
        // that is why we need to convert MPI_Request to AMPI_Request above
        PstreamGlobals::outstandingRequests_.push_back(ampiReq);
    }
}

// * * * * * * * * * * * * * * Convenience Methods * * * * * * * * * * * * * //

//- Broadcast a single int64 value.
//
//  Ensures consistent data types. Used within the following:
//  - UIPBstream::bufferIPCrecv()
//  - UOPBstream::bufferIPCsend()
//  - UOPBstream::send(Foam::zero, ...)

inline bool broadcast_int64(int64_t& value, int comm)
{
    return
    (
        MPI_SUCCESS
     == MPI_Bcast(&value, 1, MPI_INT64_T, 0, MPICommunicators_[comm])
    );
}


// * * * * * * * * * * * * * * codi: Request Helpers * * * * * * * * * * * * * * * //
// these helper functions allow us to operate based on AMPI_Request with minimal
// change to the original MPI codes

//- codi: Wait for request at index (handles both MPI and AMPI automatically)
inline void waitRequestAtIndex(const label i)
{
    // All requests (both MPI and AMPI) are now stored in outstandingRequests_
    // MPI requests use only the .request member, AMPI requests use the full structure
    if (i >= 0 && i < outstandingRequests_.size())
    {
        auto& request = outstandingRequests_[i];
        if (request.request != MPI_REQUEST_NULL)
        {
            // Check if this is a real AMPI request or just a wrapped MPI request
            // Real AMPI requests have handle != nullptr
            if (request.handle != nullptr)
            {
                // Real AMPI request: use AMPI_Wait
                if (AMPI_Wait(&request, AMPI_STATUS_IGNORE))
                {
                    FatalErrorInFunction
                        << "AMPI_Wait returned with error for index " << i
                        << Foam::abort(FatalError);
                }
            }
            else
            {
                // Wrapped MPI request: use standard MPI_Wait
                if (MPI_Wait(&request.request, MPI_STATUS_IGNORE))
                {
                    FatalErrorInFunction
                        << "MPI_Wait returned with error for index " << i
                        << Foam::abort(FatalError);
                }
            }
        }
    }
}

//- codi: Test if request at index is finished (handles both MPI and AMPI)
inline bool testRequestAtIndex(const label i)
{
    // All requests (both MPI and AMPI) are stored in outstandingRequests_
    if (i >= 0 && i < outstandingRequests_.size())
    {
        auto& request = outstandingRequests_[i];
        if (request.request == MPI_REQUEST_NULL)
        {
            return true;
        }

        int flag = 0;

        // Check if this is a real AMPI request or just a wrapped MPI request
        if (request.handle != nullptr)
        {
            // Real AMPI request: use AMPI_Test
            AMPI_Test(&request, &flag, AMPI_STATUS_IGNORE);
        }
        else
        {
            // Wrapped MPI request: use standard MPI_Test
            MPI_Test(&request.request, &flag, MPI_STATUS_IGNORE);
        }

        return flag != 0;
    }

    return true;  // Out of range or invalid index
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PstreamGlobals
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
