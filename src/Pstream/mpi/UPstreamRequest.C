/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2023-2025 OpenCFD Ltd.
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

#include "UPstreamWrapping.H"
#include "PstreamGlobals.H"
#include "profilingPstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPstream::Request::Request() noexcept
:
    UPstream::Request(MPI_REQUEST_NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UPstream::Request::good() const noexcept
{
    return MPI_REQUEST_NULL != PstreamUtils::Cast::to_mpi(*this);
}


void Foam::UPstream::Request::reset() noexcept
{
    *this = UPstream::Request(MPI_REQUEST_NULL);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// Foam::UPstream::Request
// Foam::UPstream::Request::lookup(const label req)
// {
//     if (req < 0 || req >= PstreamGlobals::outstandingRequests_.size())
//     {
//         WarningInFunction
//             << "Illegal request " << req << nl
//             << "Should be within range [0,"
//             << PstreamGlobals::outstandingRequests_.size()
//             << ')' << endl;
//
//         return UPstream::Request(MPI_REQUEST_NULL);
//     }
//
//     return UPstream::Request(PstreamGlobals::outstandingRequests_[req]);
// }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::UPstream::nRequests() noexcept
{
    return PstreamGlobals::outstandingRequests_.size();
}


void Foam::UPstream::resetRequests(const label n)
{
    if (n >= 0 && n < PstreamGlobals::outstandingRequests_.size())
    {
        PstreamGlobals::outstandingRequests_.resize(n);
    }
}


void Foam::UPstream::addRequest(UPstream::Request& req)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return;
    }

    {
        MPI_Request request = PstreamUtils::Cast::to_mpi(req);
        if (MPI_REQUEST_NULL != request)
        {
            // codi: Wrap MPI_Request in AMPI_Request structure for unified storage
            AMPI_Request ampiReq;
            ampiReq.request = request;
            ampiReq.handle = nullptr;
            ampiReq.func = nullptr;
            ampiReq.start = nullptr;
            ampiReq.end = nullptr;
            ampiReq.isActive = false;
            ampiReq.reverseData = nullptr;
            ampiReq.deleteDataFunc = nullptr;

            PstreamGlobals::outstandingRequests_.push_back(ampiReq);
        }
    }

    // Invalidate parameter
    req = UPstream::Request(MPI_REQUEST_NULL);
}


void Foam::UPstream::cancelRequest(const label i)
{
    // No-op for non-parallel, or out-of-range (eg, placeholder indices)
    if
    (
        !UPstream::parRun()
     || i < 0
     || i >= PstreamGlobals::outstandingRequests_.size()
    )
    {
        return;
    }

    {
        // codi: we use AMPI_Request
        auto& request = PstreamGlobals::outstandingRequests_[i];
        if (MPI_REQUEST_NULL != request.request)  // Active handle is mandatory
        {
            MPI_Cancel(&request.request);
            MPI_Request_free(&request.request);  //<- Sets to MPI_REQUEST_NULL
        }
    }
}


void Foam::UPstream::cancelRequest(UPstream::Request& req)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return;
    }

    {
        MPI_Request request = PstreamUtils::Cast::to_mpi(req);
        if (MPI_REQUEST_NULL != request)  // Active handle is mandatory
        {
            MPI_Cancel(&request);
            MPI_Request_free(&request);
        }
        req = UPstream::Request(MPI_REQUEST_NULL);  // Now inactive
    }
}


void Foam::UPstream::cancelRequests(UList<UPstream::Request>& requests)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return;
    }

    for (auto& req : requests)
    {
        MPI_Request request = PstreamUtils::Cast::to_mpi(req);
        if (MPI_REQUEST_NULL != request)  // Active handle is mandatory
        {
            MPI_Cancel(&request);
            MPI_Request_free(&request);
        }
        req = UPstream::Request(MPI_REQUEST_NULL);  // Now inactive
    }
}


void Foam::UPstream::removeRequests(const label pos, label len)
{
    // No-op for non-parallel, no pending requests or out-of-range
    if
    (
        !UPstream::parRun()
     || (pos < 0 || pos >= PstreamGlobals::outstandingRequests_.size())
     || !len
    )
    {
        return;
    }

    label count = (PstreamGlobals::outstandingRequests_.size() - pos);

    // Apply range-checking on slice with (len < 0) behaving like npos
    // (ie, the rest of the list)
    if (len >= 0 && len < count)
    {
        // A non-trailing slice
        count = len;
    }
    // Have count >= 1

    const labelRange range(pos, count);

    for (const label i : range)
    {
        // codi: we use AMPI_Request
        auto& request = PstreamGlobals::outstandingRequests_[i];
        if (MPI_REQUEST_NULL != request.request)  // Active handle is mandatory
        {
            MPI_Cancel(&request.request);
            MPI_Request_free(&request.request);  //<- Sets to MPI_REQUEST_NULL
        }
    }

    // Remove from list of outstanding requests and move down
    PstreamGlobals::outstandingRequests_.remove(range);
}


void Foam::UPstream::freeRequest(UPstream::Request& req)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return;
    }

    {
        MPI_Request request = PstreamUtils::Cast::to_mpi(req);
        if (MPI_REQUEST_NULL != request)  // Active handle is mandatory
        {
            // if (cancel)
            // {
            //     MPI_Cancel(&request);
            // }
            MPI_Request_free(&request);
        }
        req = UPstream::Request(MPI_REQUEST_NULL);  // Now inactive
    }
}


void Foam::UPstream::freeRequests(UList<UPstream::Request>& requests)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return;
    }

    for (auto& req : requests)
    {
        MPI_Request request = PstreamUtils::Cast::to_mpi(req);
        if (MPI_REQUEST_NULL != request)  // Active handle is mandatory
        {
            // if (cancel)
            // {
            //     MPI_Cancel(&request);
            // }
            MPI_Request_free(&request);
        }
        req = UPstream::Request(MPI_REQUEST_NULL);  // Now inactive
    }
}


void Foam::UPstream::waitRequests(const label pos, label len)
{
    // No-op for non-parallel, no pending requests or out-of-range
    if
    (
        !UPstream::parRun()
     || (pos < 0 || pos >= PstreamGlobals::outstandingRequests_.size())
     || !len
    )
    {
        return;
    }

    label count = (PstreamGlobals::outstandingRequests_.size() - pos);
    bool trim = true;  // Can trim the trailing part of the list

    // Apply range-checking on slice with (len < 0) behaving like npos
    // (ie, the rest of the list)
    if (len >= 0 && len < count)
    {
        // A non-trailing slice
        count = len;
        trim = false;
    }
    // Have count >= 1

    auto* waitRequests = (PstreamGlobals::outstandingRequests_.data() + pos);

    if (UPstream::debug)
    {
        Perr<< "UPstream::waitRequests : starting wait for "
            << count << " requests starting at " << pos << endl;
    }

    profilingPstream::beginTiming();

    if (count == 1)
    {
        // On success: sets request to AMPI_REQUEST_NULL
        if (AMPI_Wait(waitRequests, MPI_STATUS_IGNORE))
        {
            FatalErrorInFunction
                << "AMPI_Wait returned with error"
                << Foam::abort(FatalError);
        }
    }
    else if (count > 1)
    {
        // codi: On success: sets each request to AMPI_REQUEST_NULL
        // NOTE: we need to use AMPI_Wait
        if (AMPI_Waitall(count, waitRequests, MPI_STATUSES_IGNORE))
        {
            FatalErrorInFunction
                << "AMPI_Waitall returned with error"
                << Foam::abort(FatalError);
        }
    }

    profilingPstream::addWaitTime();

    if (trim)
    {
        // Trim the length of outstanding requests
        PstreamGlobals::outstandingRequests_.resize(pos);
    }

    if (UPstream::debug)
    {
        Perr<< "UPstream::waitRequests : finished wait." << endl;
    }
}


void Foam::UPstream::waitRequests(UList<UPstream::Request>& requests)
{
    // No-op for non-parallel or no pending requests
    if (!UPstream::parRun() || requests.empty())
    {
        return;
    }

    // Looks ugly but is legitimate since UPstream::Request is an intptr_t,
    // which is always large enough to hold an MPI_Request (int or pointer)

    label count = 0;
    auto* waitRequests = reinterpret_cast<MPI_Request*>(requests.data());

    for (auto& req : requests)
    {
        MPI_Request request = PstreamUtils::Cast::to_mpi(req);

        if (MPI_REQUEST_NULL != request)  // Apply some prefiltering
        {
            waitRequests[count] = request;
            ++count;
        }
    }

    if (!count)
    {
        // No active request handles
        return;
    }

    profilingPstream::beginTiming();

    // On success: sets each request to MPI_REQUEST_NULL
    if (MPI_Waitall(count, waitRequests, MPI_STATUSES_IGNORE))
    {
        FatalErrorInFunction
            << "MPI_Waitall returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    // Everything handled, reset all to MPI_REQUEST_NULL
    requests = UPstream::Request(MPI_REQUEST_NULL);
}


bool Foam::UPstream::waitAnyRequest(const label pos, label len)
{
    // No-op for non-parallel, no pending requests or out-of-range
    if
    (
        !UPstream::parRun()
     || (pos < 0 || pos >= PstreamGlobals::outstandingRequests_.size())
     || !len
    )
    {
        return false;
    }

    label count = (PstreamGlobals::outstandingRequests_.size() - pos);

    // Apply range-checking on slice with (len < 0) behaving like npos
    // (ie, the rest of the list)
    if (len >= 0 && len < count)
    {
        // A non-trailing slice
        count = len;
    }
    // Have count >= 1

    auto* waitRequests = (PstreamGlobals::outstandingRequests_.data() + pos);

    if (UPstream::debug)
    {
        Perr<< "UPstream::waitAnyRequest : starting wait for any of "
            << count << " requests starting at " << pos << endl;
    }

    profilingPstream::beginTiming();

    // On success: sets request to AMPI_REQUEST_NULL
    // codi: we need to use AMPI_Wait!
    int index = MPI_UNDEFINED;
    if (AMPI_Waitany(count, waitRequests, &index, MPI_STATUS_IGNORE))
    {
        FatalErrorInFunction
            << "AMPI_Waitany returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    if (index == MPI_UNDEFINED)
    {
        // No active request handles
        return false;
    }

    return true;
}


bool Foam::UPstream::waitSomeRequests
(
    const label pos,
    label len,
    DynamicList<int>* indices
)
{
    // No-op for non-parallel, no pending requests or out-of-range
    if
    (
        !UPstream::parRun()
     || (pos < 0 || pos >= PstreamGlobals::outstandingRequests_.size())
     || !len
    )
    {
        if (indices) indices->clear();
        return false;
    }

    label count = (PstreamGlobals::outstandingRequests_.size() - pos);

    // Apply range-checking on slice with (len < 0) behaving like npos
    // (ie, the rest of the list)
    if (len >= 0 && len < count)
    {
        // A non-trailing slice
        count = len;
    }
    // Have count >= 1

    auto* waitRequests = (PstreamGlobals::outstandingRequests_.data() + pos);

    if (UPstream::debug)
    {
        Perr<< "UPstream:waitSomeRequest : starting wait for some of "
            << count << " requests starting at " << pos << endl;
    }


    // Local temporary storage, or return via calling parameter
    List<int> tmpIndices;
    if (indices)
    {
        indices->resize_nocopy(count);
    }
    else
    {
        tmpIndices.resize(count);
    }

    profilingPstream::beginTiming();

    // On success: sets non-blocking requests to AMPI_REQUEST_NULL
    // codi: We need to use AMPI_Wait
    int outcount = 0;
    if
    (
        AMPI_Waitsome
        (
            count,
            waitRequests,
           &outcount,
            (indices ? indices->data() : tmpIndices.data()),
            MPI_STATUSES_IGNORE
        )
    )
    {
        FatalErrorInFunction
            << "AMPI_Waitsome returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    if (outcount == MPI_UNDEFINED || outcount < 1)
    {
        // No active request handles
        if (indices) indices->clear();
        return false;
    }

    if (indices)
    {
        indices->resize(outcount);
    }

    return true;
}


bool Foam::UPstream::waitSomeRequests
(
    UList<UPstream::Request>& requests,
    DynamicList<int>* indices
)
{
    // No-op for non-parallel or no pending requests
    if (!UPstream::parRun() || requests.empty())
    {
        if (indices) indices->clear();
        return false;
    }

    // Looks ugly but is legitimate since UPstream::Request is an intptr_t,
    // which is always large enough to hold an MPI_Request (int or pointer)

    label count = 0;
    auto* waitRequests = reinterpret_cast<MPI_Request*>(requests.data());

    for (auto& req : requests)
    {
        waitRequests[count] = PstreamUtils::Cast::to_mpi(req);
        ++count;
    }

    // Local temporary storage, or return via calling parameter
    List<int> tmpIndices;
    if (indices)
    {
        indices->resize_nocopy(count);
    }
    else
    {
        tmpIndices.resize(count);
    }

    if (UPstream::debug)
    {
        Perr<< "UPstream:waitSomeRequest : starting wait for some of "
            << requests.size() << " requests" << endl;
    }

    profilingPstream::beginTiming();

    // On success: sets non-blocking requests to MPI_REQUEST_NULL
    int outcount = 0;
    if
    (
        MPI_Waitsome
        (
            count,
            waitRequests,
           &outcount,
            (indices ? indices->data() : tmpIndices.data()),
            MPI_STATUSES_IGNORE
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Waitsome returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    if (outcount == MPI_UNDEFINED || outcount < 1)
    {
        // No active request handles
        if (indices) indices->clear();

        // Everything handled or inactive, reset all to MPI_REQUEST_NULL
        requests = UPstream::Request(MPI_REQUEST_NULL);
        return false;
    }

    if (indices)
    {
        indices->resize(outcount);
    }

    // Transcribe MPI_Request back into UPstream::Request
    // - do in reverse order - see note in finishedRequests()
    {
        for (label i = requests.size()-1; i >= 0; --i)
        {
            requests[i] = UPstream::Request(waitRequests[i]);
        }
    }

    return true;
}


Foam::label Foam::UPstream::waitAnyRequest(UList<UPstream::Request>& requests)
{
    // No-op for non-parallel or no pending requests
    if (!UPstream::parRun() || requests.empty())
    {
        return -1;
    }

    // Looks ugly but is legitimate since UPstream::Request is an intptr_t,
    // which is always large enough to hold an MPI_Request (int or pointer)

    label count = 0;
    auto* waitRequests = reinterpret_cast<MPI_Request*>(requests.data());

    // Transcribe UPstream::Request into MPI_Request
    // - do not change locations within the list since these are relevant
    //   for the return index.
    for (auto& req : requests)
    {
        waitRequests[count] = PstreamUtils::Cast::to_mpi(req);
        ++count;
    }

    profilingPstream::beginTiming();

    // On success: sets request to MPI_REQUEST_NULL
    int index = MPI_UNDEFINED;
    if (MPI_Waitany(count, waitRequests, &index, MPI_STATUS_IGNORE))
    {
        FatalErrorInFunction
            << "MPI_Waitany returned with error"
            << Foam::abort(FatalError);
    }

    profilingPstream::addWaitTime();

    if (index == MPI_UNDEFINED)
    {
        index = -1;  // No outstanding requests
    }

    // Transcribe MPI_Request back into UPstream::Request
    // - do in reverse order - see note in finishedRequests()
    {
        for (label i = count-1; i >= 0; --i)
        {
            requests[i] = UPstream::Request(waitRequests[i]);
        }

        // Trailing portion
        for (label i = count; i < requests.size(); ++i)
        {
            requests[i] = UPstream::Request(MPI_REQUEST_NULL);
        }
    }

    return index;
}


// FUTURE?
//
/// void Foam::UPstream::waitRequests
/// (
///     UPstream::Request& req0,
///     UPstream::Request& req1
/// )
/// {
///     // No-op for non-parallel
///     if (!UPstream::parRun())
///     {
///         return;
///     }
///
///     int count = 0;
///     MPI_Request waitRequests[2];
///
///     waitRequests[count] = PstreamUtils::Cast::to_mpi(req0);
///     if (MPI_REQUEST_NULL != waitRequests[count])
///     {
///         ++count;
///     }
///
///     waitRequests[count] = PstreamUtils::Cast::to_mpi(req1);
///     if (MPI_REQUEST_NULL != waitRequests[count])
///     {
///         ++count;
///     }
///
///     // Flag in advance as being handled
///     req0 = UPstream::Request(MPI_REQUEST_NULL);
///     req1 = UPstream::Request(MPI_REQUEST_NULL);
///
///     if (!count)
///     {
///         return;
///     }
///
///     profilingPstream::beginTiming();
///
///     // On success: sets each request to MPI_REQUEST_NULL
///     if (MPI_Waitall(count, waitRequests, MPI_STATUSES_IGNORE))
///     {
///         FatalErrorInFunction
///             << "MPI_Waitall returned with error"
///             << Foam::abort(FatalError);
///     }
///
///     profilingPstream::addWaitTime();
/// }


void Foam::UPstream::waitRequest(const label i)
{
    // No-op for non-parallel or invalid index
    // codi: we have to remove the || condition for the size check because
    // we save both MPI and AMPI requests
    if (!UPstream::parRun() || i < 0)
    {
        return;
    }

    if (UPstream::debug)
    {
        Perr<< "UPstream::waitRequest : starting wait for request:"
            << i << endl;
    }

    profilingPstream::beginTiming();

    // codi: Use helper that handles both MPI and AMPI automatically
    PstreamGlobals::waitRequestAtIndex(i);

    profilingPstream::addWaitTime();

    if (UPstream::debug)
    {
        Perr<< "UPstream::waitRequest : finished wait for request:"
            << i << endl;
    }
}


void Foam::UPstream::waitRequest(UPstream::Request& req)
{
    // codi: we need to do major change to the waitRequest code logic for AMPI_Request

    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return;
    }

    // codi: get the request value
    std::intptr_t value = req.value();

    // codi: No-op for null request
    if (value == 0)
    {
        return;
    }

    profilingPstream::beginTiming();

    // codi: Check if this is an AD request (negative index) or MPI request (positive)
    if (value < 0)
    {
        // AD request: stored as negative index
        // Decode: -1 → index 0, -2 → index 1, etc.
        label index = -(value + 1);

        if (index >= 0 && index < PstreamGlobals::outstandingRequests_.size())
        {
            AMPI_Request& request = PstreamGlobals::outstandingRequests_[index];

            if (AMPI_Wait(&request, MPI_STATUS_IGNORE))
            {
                FatalErrorInFunction
                    << "AMPI_Wait returned with error"
                    << Foam::abort(FatalError);
            }

            // Mark as completed (set to REQUEST_NULL)
            request.request = MPI_REQUEST_NULL;
        }
        else
        {
            FatalErrorInFunction
                << "Invalid AD request index: " << index
                << " (size=" << PstreamGlobals::outstandingRequests_.size() << ")"
                << Foam::abort(FatalError);
        }
    }
    else
    {
        // Standard MPI request: stored as positive value
        MPI_Request request = PstreamUtils::Cast::to_mpi(req);

        if (MPI_Wait(&request, MPI_STATUS_IGNORE))
        {
            FatalErrorInFunction
                << "MPI_Wait returned with error"
                << Foam::abort(FatalError);
        }
    }

    profilingPstream::addWaitTime();

    req = UPstream::Request(MPI_REQUEST_NULL);  // Now inactive
}


bool Foam::UPstream::finishedRequest(const label i)
{
    // codi: No-op for non-parallel or invalid index
    // we have to remove the || condition for the size check because
    // we save both MPI and AMPI requests
    if (!UPstream::parRun() || i < 0)
    {
        return true;
    }

    if (UPstream::debug)
    {
        Perr<< "UPstream::finishedRequest : check request:"
            << i << endl;
    }

    // codi: Use helper that handles both MPI and AMPI automatically
    return PstreamGlobals::testRequestAtIndex(i);
}


bool Foam::UPstream::finishedRequest(UPstream::Request& req)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        return true;
    }

    // codi:
    std::intptr_t value = req.value();

    // codi: Fast-path (no-op) for null request
    if (value == 0)
    {
        return true;
    }

    int flag = 0;

    // codi: Check if this is an AD request (negative index) or MPI request (positive)
    if (value < 0)
    {
        // AD request: stored as negative index
        // Decode: -1 → index 0, -2 → index 1, etc.
        label index = -(value + 1);

        if (index >= 0 && index < PstreamGlobals::outstandingRequests_.size())
        {
            AMPI_Request& request = PstreamGlobals::outstandingRequests_[index];

            // Fast-path for already completed request
            if (MPI_REQUEST_NULL == request.request)
            {
                return true;
            }

            AMPI_Test(&request, &flag, MPI_STATUS_IGNORE);

            if (flag)
            {
                // Success: mark as completed
                request.request = MPI_REQUEST_NULL;
                req = UPstream::Request(MPI_REQUEST_NULL);
            }
        }
        else
        {
            FatalErrorInFunction
                << "Invalid AD request index: " << index
                << " (size=" << PstreamGlobals::outstandingRequests_.size() << ")"
                << Foam::abort(FatalError);
        }
    }
    else
    {
        // Standard MPI request: stored as positive value
        MPI_Request request = PstreamUtils::Cast::to_mpi(req);

        MPI_Test(&request, &flag, MPI_STATUS_IGNORE);

        if (flag)
        {
            // Success: now inactive
            req = UPstream::Request(MPI_REQUEST_NULL);
        }
    }

    return flag != 0;
}


bool Foam::UPstream::finishedRequests(const label pos, label len)
{
    // No-op for non-parallel, or out-of-range (eg, placeholder indices)
    if
    (
        !UPstream::parRun()
     || (pos < 0 || pos >= PstreamGlobals::outstandingRequests_.size())
     || !len
    )
    {
        return true;
    }

    label count = (PstreamGlobals::outstandingRequests_.size() - pos);

    // Apply range-checking on slice with (len < 0) behaving like npos
    // (ie, the rest of the list)
    if (len >= 0 && len < count)
    {
        // A non-trailing slice
        count = len;
    }
    // Have count >= 1

    if (UPstream::debug)
    {
        Perr<< "UPstream::finishedRequests : check " << count
            << " requests starting at " << pos << endl;
    }

    auto* waitRequests = (PstreamGlobals::outstandingRequests_.data() + pos);

    int flag = 1;

    if (count == 1)
    {
        // Fast-path (no-op) for single null request
        // codi:
        if (MPI_REQUEST_NULL == waitRequests->request)
        {
            return true;
        }

        // On success: sets request to AMPI_REQUEST_NULL
        // codi:
        AMPI_Test(waitRequests, &flag, MPI_STATUS_IGNORE);
    }
    else if (count > 1)
    {
        // On success: sets each request to AMPI_REQUEST_NULL
        // On failure: no request is modified
        // codi:
        AMPI_Testall(count, waitRequests, &flag, MPI_STATUSES_IGNORE);
    }

    return flag != 0;
}


bool Foam::UPstream::finishedRequests(UList<UPstream::Request>& requests)
{
    // No-op for non-parallel or no pending requests
    if (!UPstream::parRun() || requests.empty())
    {
        return true;
    }

    // Looks ugly but is legitimate since UPstream::Request is an intptr_t,
    // which is always large enough to hold an MPI_Request (int or pointer)

    label count = 0;
    auto* waitRequests = reinterpret_cast<MPI_Request*>(requests.data());

    for (auto& req : requests)
    {
        MPI_Request request = PstreamUtils::Cast::to_mpi(req);

        if (MPI_REQUEST_NULL != request)  // Apply some prefiltering
        {
            waitRequests[count] = request;
            ++count;
        }
    }

    if (!count)
    {
        // No active handles
        return true;
    }

    // On success: sets each request to MPI_REQUEST_NULL
    // On failure: no request is modified
    int flag = 0;
    MPI_Testall(count, waitRequests, &flag, MPI_STATUSES_IGNORE);

    if (flag)
    {
        // Success: reset all requests to MPI_REQUEST_NULL
        requests = UPstream::Request(MPI_REQUEST_NULL);
    }
    else
    {
        // Not all done. Recover wrapped representation but in reverse order
        // since sizeof(MPI_Request) can be smaller than
        // sizeof(UPstream::Request::value_type)
        // eg, mpich has MPI_Request as 'int'
        //
        // This is uglier that we'd like, but much better than allocating
        // and freeing a scratch buffer each time we query things.

        for (label i = count-1; i >= 0; --i)
        {
            requests[i] = UPstream::Request(waitRequests[i]);
        }

        // Trailing portion
        for (label i = count; i < requests.size(); ++i)
        {
            requests[i] = UPstream::Request(MPI_REQUEST_NULL);
        }
    }

    return flag != 0;
}


bool Foam::UPstream::finishedRequestPair(label& req0, label& req1)
{
    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        req0 = -1;
        req1 = -1;
        return true;
    }

    // codi: Use helper functions to test each request individually
    bool finished0 = PstreamGlobals::testRequestAtIndex(req0);
    bool finished1 = PstreamGlobals::testRequestAtIndex(req1);

    // Mark finished requests as done
    if (finished0)
    {
        req0 = -1;
    }

    if (finished1)
    {
        req1 = -1;
    }

    // Return true if both are finished
    return (finished0 && finished1);
}


void Foam::UPstream::waitRequestPair(label& req0, label& req1)
{
    // No-op for non-parallel. Flag indices as 'done'
    if (!UPstream::parRun())
    {
        req0 = -1;
        req1 = -1;
        return;
    }

    // codi: Use helper functions to wait for each request individually
    PstreamGlobals::waitRequestAtIndex(req0);
    PstreamGlobals::waitRequestAtIndex(req1);

    // Flag indices as 'done'
    req0 = -1;
    req1 = -1;
}


// ************************************************************************* //
