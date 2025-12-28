

#include "fvCFD.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"

    codi::RealReverse::Tape& tape = codi::RealReverse::getTape();

    // Track test failures
    bool allTestsPassed = true;

    // ************************************************************************
    // Test 1: AD scalar reduce with sumOp
    //
    // Description: Tests AMPI_Allreduce with sumOp on AD scalar types
    // Each rank computes val^2 where val = rank + 3, then reduces sum
    // Expectation: Derivatives propagate correctly - df/dval = 2*val for each rank
    // ************************************************************************
    {
        tape.reset();
        tape.setActive();

        scalar val = Pstream::myProcNo() + 3;
        tape.registerInput(val);

        scalar val2 = val * val;
        reduce(val2, sumOp<scalar>());

        tape.registerOutput(val2);
        tape.setPassive();

        // CRITICAL: Only master rank sets gradient to avoid MediPack summing across ranks
        if (Pstream::master())
        {
            val2.setGradient(1.0);
        }

        // Ensure all ranks are synchronized before tape evaluation
        tape.evaluate();

        scalar expected_deriv = 2.0 * val.getValue();
        scalar computed_deriv = val.getGradient();
        bool passed = mag(computed_deriv - expected_deriv) < 1e-10;
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 1 (AD scalar sumOp): fVal=" << val2
             << " df/dval=" << computed_deriv
             << " (expected=" << expected_deriv << ")"
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // ************************************************************************
    // Test 2: AD scalar reduce with maxOp
    //
    // Description: Tests AMPI_Allreduce with maxOp on AD scalar types
    // Each rank computes val^2 where val = rank + 3, then reduces max
    // Expectation: Only the rank with max value gets derivative = 2*val
    //              (highest rank nProcs-1), all others get 0
    // ************************************************************************
    {
        tape.reset();
        tape.setActive();

        scalar val = Pstream::myProcNo() + 3;
        tape.registerInput(val);

        scalar val2 = val * val;
        reduce(val2, maxOp<scalar>());

        tape.registerOutput(val2);
        tape.setPassive();

        // CRITICAL: Only master rank sets gradient to avoid MediPack summing across ranks
        if (Pstream::master())
        {
            val2.setGradient(1.0);
        }

        // Ensure all ranks are synchronized before tape evaluation
        tape.evaluate();

        int maxRank = Pstream::nProcs() - 1;
        scalar expected_deriv = (Pstream::myProcNo() == maxRank) ? 2.0 * val.getValue() : 0.0;
        scalar computed_deriv = val.getGradient();
        bool passed = mag(computed_deriv - expected_deriv) < 1e-10;
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 2 (AD scalar maxOp): fVal=" << val2
             << " df/dval=" << computed_deriv
             << " (expected=" << expected_deriv << ")"
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // ************************************************************************
    // Test 3: AD scalar reduce with minOp
    //
    // Description: Tests AMPI_Allreduce with minOp on AD scalar types
    // Each rank computes val^2 where val = rank + 3, then reduces min
    // Expectation: Only rank 0 (min value) gets derivative = 2*val,
    //              all others get 0
    // ************************************************************************
    {
        tape.reset();
        tape.setActive();

        scalar val = Pstream::myProcNo() + 3;
        tape.registerInput(val);

        scalar val2 = val * val;
        reduce(val2, minOp<scalar>());

        tape.registerOutput(val2);
        tape.setPassive();

        // CRITICAL: Only master rank sets gradient to avoid MediPack summing across ranks
        if (Pstream::master())
        {
            val2.setGradient(1.0);
        }
        tape.evaluate();

        auto expected_deriv = (Pstream::myProcNo() == 0) ? 2.0 * val : 0.0 * val;
        auto computed_deriv = val.getGradient();
        bool passed = std::abs(computed_deriv - expected_deriv) < 1e-10;
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 3 (AD scalar minOp): fVal=" << val2
             << " df/dval=" << computed_deriv
             << " (expected=" << expected_deriv << ")"
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // ************************************************************************
    // Test 4: AD vector reduce with sumOp
    //
    // Description: Tests AMPI_Allreduce with sumOp on AD vector types
    // Each rank has vector v = (rank+1, rank+2, rank+3), computes v^2 componentwise
    // Expectation: Each component derivative = 2*component_value for each rank
    // ************************************************************************
    {
        tape.reset();
        tape.setActive();

        vector v
        (
            Pstream::myProcNo() + 1,
            Pstream::myProcNo() + 2,
            Pstream::myProcNo() + 3
        );

        tape.registerInput(v.x());
        tape.registerInput(v.y());
        tape.registerInput(v.z());

        vector v2(v.x() * v.x(), v.y() * v.y(), v.z() * v.z());
        reduce(v2, sumOp<vector>());

        tape.registerOutput(v2.x());
        tape.registerOutput(v2.y());
        tape.registerOutput(v2.z());
        tape.setPassive();

        // CRITICAL: Only master rank sets gradients to avoid MediPack summing across ranks
        if (Pstream::master())
        {
            v2.x().setGradient(1.0);
            v2.y().setGradient(1.0);
            v2.z().setGradient(1.0);
        }

        tape.evaluate();

        vector expected_deriv
        (
            2.0 * v.x().getValue(),
            2.0 * v.y().getValue(),
            2.0 * v.z().getValue()
        );

        vector computed_deriv
        (
            v.x().getGradient(),
            v.y().getGradient(),
            v.z().getGradient()
        );

        bool passed = mag(computed_deriv - expected_deriv) < 1e-10;
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 4 (AD vector sumOp): fVal=" << v2
             << " df/dv=" << computed_deriv
             << " (expected=" << expected_deriv << ")"
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // ************************************************************************
    // Test 5: AD vector reduce with maxOp
    //
    // Description: Tests AMPI_Allreduce with maxOp on AD vector types
    // Vectors have increasing magnitude: rank 0: (1,2,3), rank 1: (2,3,4), etc.
    // Expectation: Only highest rank gets derivatives (has max magnitude)
    // ************************************************************************
    {
        tape.reset();
        tape.setActive();

        vector v
        (
            Pstream::myProcNo() + 1,
            Pstream::myProcNo() + 2,
            Pstream::myProcNo() + 3
        );

        tape.registerInput(v.x());
        tape.registerInput(v.y());
        tape.registerInput(v.z());

        vector v2(v.x() * v.x(), v.y() * v.y(), v.z() * v.z());
        reduce(v2, maxOp<vector>());

        tape.registerOutput(v2.x());
        tape.registerOutput(v2.y());
        tape.registerOutput(v2.z());
        tape.setPassive();

        // CRITICAL: Only master rank sets gradients
        if (Pstream::master())
        {
            v2.x().setGradient(1.0);
            v2.y().setGradient(1.0);
            v2.z().setGradient(1.0);
        }

        tape.evaluate();

        int maxRank = Pstream::nProcs() - 1;
        vector expected_deriv = vector::zero;
        if (Pstream::myProcNo() == maxRank)
        {
            expected_deriv.x() = 2.0 * v.x().getValue();
            expected_deriv.y() = 2.0 * v.y().getValue();
            expected_deriv.z() = 2.0 * v.z().getValue();
        }

        vector computed_deriv
        (
            v.x().getGradient(),
            v.y().getGradient(),
            v.z().getGradient()
        );

        bool passed = mag(computed_deriv - expected_deriv) < 1e-10;
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 5 (AD vector maxOp): fVal=" << v2
             << " df/dv=" << computed_deriv
             << " (expected=" << expected_deriv << ")"
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // ************************************************************************
    // Test 6: AD vector reduce with minOp
    //
    // Description: Tests AMPI_Allreduce with minOp on AD vector types
    // Vectors have increasing magnitude
    // Expectation: Only rank 0 (min magnitude) gets derivatives
    // ************************************************************************
    {
        tape.reset();
        tape.setActive();

        vector v
        (
            Pstream::myProcNo() + 1,
            Pstream::myProcNo() + 2,
            Pstream::myProcNo() + 3
        );

        tape.registerInput(v.x());
        tape.registerInput(v.y());
        tape.registerInput(v.z());

        vector v2(v.x() * v.x(), v.y() * v.y(), v.z() * v.z());
        reduce(v2, minOp<vector>());

        tape.registerOutput(v2.x());
        tape.registerOutput(v2.y());
        tape.registerOutput(v2.z());
        tape.setPassive();

        // CRITICAL: Only master rank sets gradients
        if (Pstream::master())
        {
            v2.x().setGradient(1.0);
            v2.y().setGradient(1.0);
            v2.z().setGradient(1.0);
        }

        tape.evaluate();

        vector expected_deriv = vector::zero;
        if (Pstream::myProcNo() == 0)
        {
            expected_deriv.x() = 2.0 * v.x().getValue();
            expected_deriv.y() = 2.0 * v.y().getValue();
            expected_deriv.z() = 2.0 * v.z().getValue();
        }

        vector computed_deriv
        (
            v.x().getGradient(),
            v.y().getGradient(),
            v.z().getGradient()
        );

        bool passed = mag(computed_deriv - expected_deriv) < 1e-10;
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 6 (AD vector minOp): fVal=" << v2
             << " df/dv=" << computed_deriv
             << " (expected=" << expected_deriv << ")"
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // ************************************************************************
    // Test 7: Non-AD label reduce with sumOp
    //
    // Description: Tests standard MPI_Allreduce (non-AD) with sumOp on label (integer)
    // Each rank contributes rank+1: rank 0→1, rank 1→2, rank 2→3, rank 3→4
    // Expectation: Sum = nProcs*(nProcs+1)/2 (e.g., 4 ranks: 1+2+3+4=10)
    //              No AD involved, tests standard MPI path
    // ************************************************************************
    {
        int myRank = Pstream::myProcNo();
        int nProcs = Pstream::nProcs();

        label intVal = myRank + 1;
        label originalVal = intVal;

        reduce(intVal, sumOp<label>());

        // Expected sum: 1 + 2 + 3 + ... + nProcs = nProcs*(nProcs+1)/2
        label expected_sum = nProcs * (nProcs + 1) / 2;
        bool passed = (intVal == expected_sum);
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 7 (non-AD label sumOp): myVal=" << originalVal
             << " reduced=" << intVal
             << " (expected=" << expected_sum << ")"
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // ************************************************************************
    // Test 8: AD blocking point-to-point communication with 2-element Field
    //
    // Description: Tests AMPI blocking send/recv with circular communication pattern
    // Each rank computes sendVal1^2 and sendVal2^3, sends both to next rank
    // Uses blocking buffered communication (AMPI_Send/AMPI_Recv)
    // Expectation: Derivatives propagate backwards through MPI communication
    //              Master receives from rank N-1, so rank N-1 gets derivatives
    // ************************************************************************
    {
        tape.reset();
        tape.setActive();

        int myRank = Pstream::myProcNo();
        int nProcs = Pstream::nProcs();
        int sendTo = (myRank + 1) % nProcs;
        int recvFrom = (myRank - 1 + nProcs) % nProcs;

        // Each rank creates two scalar values based on its rank
        scalar sendVal1 = myRank + 10.0;
        scalar sendVal2 = myRank + 20.0;
        tape.registerInput(sendVal1);
        tape.registerInput(sendVal2);

        // Compute functions of the inputs
        scalar sendVal1_sqr = sendVal1 * sendVal1;
        scalar sendVal2_cube = sendVal2 * sendVal2 * sendVal2;

        // Blocking send/receive using UIPstream/UOPstream
        Field<scalar> recvBuf(2);
        Field<scalar> sendBuf(2);
        sendBuf[0] = sendVal1_sqr;
        sendBuf[1] = sendVal2_cube;

        // Use blocking buffered communication
        // Send first, then receive to avoid deadlock in circular pattern
        UOPstream::write
        (
            UPstream::commsTypes::buffered,
            sendTo,
            sendBuf.cdata(),
            sendBuf.size(),
            Pstream::msgType(),
            Pstream::worldComm
        );

        UIPstream::read
        (
            UPstream::commsTypes::buffered,
            recvFrom,
            recvBuf.data(),
            recvBuf.size(),
            Pstream::msgType(),
            Pstream::worldComm
        );

        scalar recvVal1_sqr = recvBuf[0];
        scalar recvVal2_cube = recvBuf[1];

        tape.registerOutput(recvVal1_sqr);
        tape.registerOutput(recvVal2_cube);
        tape.setPassive();

        // Set gradients on the received values
        // Only master sets gradients to avoid summing
        if (Pstream::master())
        {
            recvVal1_sqr.setGradient(1.0);
            recvVal2_cube.setGradient(1.0);
        }

        tape.evaluate();

        // The gradients should propagate back to the sending rank
        // For master (rank 0), it receives from rank (nProcs-1)
        // So rank (nProcs-1)'s sendVal1 should have gradient = 2*sendVal1
        // and sendVal2 should have gradient = 3*sendVal2^2
        scalar computed_deriv1 = sendVal1.getGradient();
        scalar computed_deriv2 = sendVal2.getGradient();

        // Expected: rank (nProcs-1) should have derivatives, others 0
        scalar expected_deriv1 = 0.0;
        scalar expected_deriv2 = 0.0;
        if (myRank == nProcs - 1)
        {
            expected_deriv1 = 2.0 * sendVal1.getValue();
            expected_deriv2 = 3.0 * sendVal2.getValue() * sendVal2.getValue();
        }

        bool passed = (mag(computed_deriv1 - expected_deriv1) < 1e-10) &&
                      (mag(computed_deriv2 - expected_deriv2) < 1e-10);
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 8 (AD blocking send/recv 2-element): "
             << "sendVal1=" << sendVal1 << " sendVal2=" << sendVal2
             << " | recvVal1_sqr=" << recvVal1_sqr << " recvVal2_cube=" << recvVal2_cube
             << " | df/dsendVal1=" << computed_deriv1 << " (expected=" << expected_deriv1 << ")"
             << " df/dsendVal2=" << computed_deriv2 << " (expected=" << expected_deriv2 << ")"
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // ************************************************************************
    // Test 9: Non-AD blocking point-to-point with labelList
    //
    // Description: Tests standard MPI blocking send/recv (non-AD) with labelList
    // Each rank sends a labelList [rank*10, rank*10+1, ..., rank*10+4] to next rank
    // Uses blocking buffered communication (MPI_Send/MPI_Recv, no AD)
    // Expectation: Received data matches what previous rank sent
    //              rank 0 receives [30,31,32,33,34] from rank 3 (in 4-rank case)
    //              No AD tape involved, tests standard MPI path
    // ************************************************************************
    {
        int myRank = Pstream::myProcNo();
        int nProcs = Pstream::nProcs();
        int sendTo = (myRank + 1) % nProcs;
        int recvFrom = (myRank - 1 + nProcs) % nProcs;

        // Create a labelList with rank-specific values
        labelList sendData(5);
        for (label i = 0; i < 5; ++i)
        {
            sendData[i] = myRank * 10 + i;  // rank 0: [0,1,2,3,4], rank 1: [10,11,12,13,14], etc.
        }

        labelList recvData(5);

        // Use blocking buffered communication
        // Send first, then receive to avoid deadlock in circular pattern
        UOPstream::write
        (
            UPstream::commsTypes::buffered,
            sendTo,
            sendData.cdata(),
            sendData.size(),
            Pstream::msgType(),
            Pstream::worldComm
        );

        UIPstream::read
        (
            UPstream::commsTypes::buffered,
            recvFrom,
            recvData.data(),
            recvData.size(),
            Pstream::msgType(),
            Pstream::worldComm
        );

        // Verify received data matches what the previous rank should have sent
        label expectedBase = recvFrom * 10;
        bool passed = true;
        for (label i = 0; i < 5; ++i)
        {
            if (recvData[i] != expectedBase + i)
            {
                passed = false;
                break;
            }
        }
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 9 (non-AD blocking labelList): sent=" << sendData
             << " received=" << recvData
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // ************************************************************************
    // Test 10: AD non-blocking point-to-point communication with 2-element Field
    //
    // Description: Tests AMPI non-blocking send/recv with circular communication pattern
    // Each rank computes sendVal1^2 and sendVal2^3, sends both to next rank
    // Uses non-blocking communication (AMPI_Isend/AMPI_Irecv with request handling)
    // Expectation: Derivatives propagate backwards through MPI communication
    //              Master receives from rank N-1, so rank N-1 gets derivatives
    // ************************************************************************
    {
        tape.reset();
        tape.setActive();

        int myRank = Pstream::myProcNo();
        int nProcs = Pstream::nProcs();
        int sendTo = (myRank + 1) % nProcs;
        int recvFrom = (myRank - 1 + nProcs) % nProcs;

        // Each rank creates two scalar values based on its rank
        scalar sendVal1 = myRank + 10.0;
        scalar sendVal2 = myRank + 20.0;
        tape.registerInput(sendVal1);
        tape.registerInput(sendVal2);

        // Compute functions of the inputs
        scalar sendVal1_sqr = sendVal1 * sendVal1;
        scalar sendVal2_cube = sendVal2 * sendVal2 * sendVal2;

        // Non-blocking send/receive using UIPstream/UOPstream
        Field<scalar> recvBuf(2);
        Field<scalar> sendBuf(2);
        sendBuf[0] = sendVal1_sqr;
        sendBuf[1] = sendVal2_cube;

        // Use non-blocking communication
        UPstream::Request sendReq;
        UPstream::Request recvReq;

        // Post non-blocking receive first
        UIPstream::read
        (
            recvReq,
            recvFrom,
            recvBuf.data(),
            recvBuf.size(),
            Pstream::msgType(),
            Pstream::worldComm
        );

        // Post non-blocking send
        UOPstream::write
        (
            sendReq,
            sendTo,
            sendBuf.cdata(),
            sendBuf.size(),
            Pstream::msgType(),
            Pstream::worldComm
        );

        // Wait for both operations to complete
        UPstream::waitRequest(recvReq);
        UPstream::waitRequest(sendReq);

        scalar recvVal1_sqr = recvBuf[0];
        scalar recvVal2_cube = recvBuf[1];

        tape.registerOutput(recvVal1_sqr);
        tape.registerOutput(recvVal2_cube);
        tape.setPassive();

        // Set gradients on the received values
        // Only master sets gradients to avoid summing
        if (Pstream::master())
        {
            recvVal1_sqr.setGradient(1.0);
            recvVal2_cube.setGradient(1.0);
        }

        tape.evaluate();

        // The gradients should propagate back to the sending rank
        // For master (rank 0), it receives from rank (nProcs-1)
        // So rank (nProcs-1)'s sendVal1 should have gradient = 2*sendVal1
        // and sendVal2 should have gradient = 3*sendVal2^2
        scalar computed_deriv1 = sendVal1.getGradient();
        scalar computed_deriv2 = sendVal2.getGradient();

        // Expected: rank (nProcs-1) should have derivatives, others 0
        scalar expected_deriv1 = 0.0;
        scalar expected_deriv2 = 0.0;
        if (myRank == nProcs - 1)
        {
            expected_deriv1 = 2.0 * sendVal1.getValue();
            expected_deriv2 = 3.0 * sendVal2.getValue() * sendVal2.getValue();
        }

        bool passed = (mag(computed_deriv1 - expected_deriv1) < 1e-10) &&
                      (mag(computed_deriv2 - expected_deriv2) < 1e-10);
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 10 (AD non-blocking send/recv 2-element): "
             << "sendVal1=" << sendVal1 << " sendVal2=" << sendVal2
             << " | recvVal1_sqr=" << recvVal1_sqr << " recvVal2_cube=" << recvVal2_cube
             << " | df/dsendVal1=" << computed_deriv1 << " (expected=" << expected_deriv1 << ")"
             << " df/dsendVal2=" << computed_deriv2 << " (expected=" << expected_deriv2 << ")"
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // Test 11: Non-AD non-blocking with label-based API
    {
        tape.reset();

        int myRank = Pstream::myProcNo();
        int nProcs = Pstream::nProcs();
        int sendTo = (myRank + 1) % nProcs;
        int recvFrom = (myRank - 1 + nProcs) % nProcs;

        double sendVal = myRank + 100.0;
        double recvVal = 0.0;

        // Use label-based API: track the request indices manually
        label startIdx = UPstream::nRequests();

        UIPstream::read
        (
            UPstream::commsTypes::nonBlocking,
            recvFrom,
            reinterpret_cast<char*>(&recvVal),
            sizeof(double),
            0,  // tag
            UPstream::worldComm
        );

        UOPstream::write
        (
            UPstream::commsTypes::nonBlocking,
            sendTo,
            reinterpret_cast<const char*>(&sendVal),
            sizeof(double),
            0,  // tag
            UPstream::worldComm
        );

        // Wait using indices into the global request list
        label recvReqIdx = startIdx;      // First request pushed
        label sendReqIdx = startIdx + 1;  // Second request pushed

        UPstream::waitRequest(recvReqIdx);
        UPstream::waitRequest(sendReqIdx);

        bool passed = (recvVal == (recvFrom + 100.0));
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 11 (Non-AD waitRequest(label)): "
             << "sendVal=" << sendVal << " recvVal=" << recvVal
             << " expected=" << (recvFrom + 100.0)
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // Test 12: Non-AD non-blocking with UPstream::Request API
    {
        tape.reset();

        int myRank = Pstream::myProcNo();
        int nProcs = Pstream::nProcs();
        int sendTo = (myRank + 1) % nProcs;
        int recvFrom = (myRank - 1 + nProcs) % nProcs;

        double sendVal = myRank + 200.0;
        double recvVal = 0.0;

        // Use UPstream::Request API
        UPstream::Request recvReq;
        UPstream::Request sendReq;

        UIPstream::read
        (
            recvReq,
            recvFrom,
            reinterpret_cast<char*>(&recvVal),
            sizeof(double),
            0,  // tag
            UPstream::worldComm
        );

        UOPstream::write
        (
            sendReq,
            sendTo,
            reinterpret_cast<const char*>(&sendVal),
            sizeof(double),
            0,  // tag
            UPstream::worldComm
        );

        // Use waitRequest(UPstream::Request&) API
        UPstream::waitRequest(recvReq);
        UPstream::waitRequest(sendReq);

        bool passed = (recvVal == (recvFrom + 200.0));
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 12 (Non-AD waitRequest(Request&)): "
             << "sendVal=" << sendVal << " recvVal=" << recvVal
             << " expected=" << (recvFrom + 200.0)
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // Test 13: AD non-blocking with waitRequest(label) API
    {
        tape.reset();
        tape.setActive();

        int myRank = Pstream::myProcNo();
        int nProcs = Pstream::nProcs();
        int sendTo = (myRank + 1) % nProcs;
        int recvFrom = (myRank - 1 + nProcs) % nProcs;

        scalar sendVal = myRank + 300.0;
        tape.registerInput(sendVal);
        scalar sendVal_sqr = sendVal * sendVal;

        Field<scalar> recvBuf(1);
        Field<scalar> sendBuf(1);
        sendBuf[0] = sendVal_sqr;

        // Use label-based API: track indices manually
        label startIdx = UPstream::nRequests();

        UIPstream::read
        (
            UPstream::commsTypes::nonBlocking,
            recvFrom,
            recvBuf.data(),
            recvBuf.size(),
            0,  // tag
            UPstream::worldComm
        );

        UOPstream::write
        (
            UPstream::commsTypes::nonBlocking,
            sendTo,
            sendBuf.cdata(),
            sendBuf.size(),
            0,  // tag
            UPstream::worldComm
        );

        // Wait using calculated indices
        label recvReqIdx = startIdx;      // First request
        label sendReqIdx = startIdx + 1;  // Second request

        UPstream::waitRequest(recvReqIdx);
        UPstream::waitRequest(sendReqIdx);

        scalar recvVal_sqr = recvBuf[0];
        tape.registerOutput(recvVal_sqr);
        tape.setPassive();

        // Only master sets gradient to avoid MediPack summing across ranks
        if (Pstream::master())
        {
            recvVal_sqr.setGradient(1.0);
        }

        tape.evaluate();

        scalar computed_deriv = sendVal.getGradient();
        scalar expected_deriv = 0.0;
        if (myRank == nProcs - 1)
        {
            expected_deriv = 2.0 * sendVal.getValue();
        }

        bool passed = (mag(computed_deriv - expected_deriv) < 1e-10);
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 13 (AD waitRequest(label)): "
             << "sendVal=" << sendVal << " recvVal_sqr=" << recvVal_sqr
             << " df/dsendVal=" << computed_deriv << " expected=" << expected_deriv
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // Test 14: AD non-blocking with finishedRequest(label) API
    {
        tape.reset();
        tape.setActive();

        int myRank = Pstream::myProcNo();
        int nProcs = Pstream::nProcs();
        int sendTo = (myRank + 1) % nProcs;
        int recvFrom = (myRank - 1 + nProcs) % nProcs;

        scalar sendVal = myRank + 400.0;
        tape.registerInput(sendVal);
        scalar sendVal_cube = sendVal * sendVal * sendVal;

        Field<scalar> recvBuf(1);
        Field<scalar> sendBuf(1);
        sendBuf[0] = sendVal_cube;

        // Use label-based API: track indices manually
        label startIdx = UPstream::nRequests();

        UIPstream::read
        (
            UPstream::commsTypes::nonBlocking,
            recvFrom,
            recvBuf.data(),
            recvBuf.size(),
            0,
            UPstream::worldComm
        );

        UOPstream::write
        (
            UPstream::commsTypes::nonBlocking,
            sendTo,
            sendBuf.cdata(),
            sendBuf.size(),
            0,
            UPstream::worldComm
        );

        // Calculate indices
        label recvReqIdx = startIdx;
        label sendReqIdx = startIdx + 1;

        // Poll with finishedRequest(label) until both complete
        while (!UPstream::finishedRequest(recvReqIdx) ||
               !UPstream::finishedRequest(sendReqIdx))
        {
            // Busy wait
        }

        scalar recvVal_cube = recvBuf[0];
        tape.registerOutput(recvVal_cube);
        tape.setPassive();

        // Only master sets gradient to avoid MediPack summing across ranks
        if (Pstream::master())
        {
            recvVal_cube.setGradient(1.0);
        }

        tape.evaluate();

        scalar computed_deriv = sendVal.getGradient();
        scalar expected_deriv = 0.0;
        if (myRank == nProcs - 1)
        {
            expected_deriv = 3.0 * sendVal.getValue() * sendVal.getValue();
        }

        bool passed = (mag(computed_deriv - expected_deriv) < 1e-10);
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 14 (AD finishedRequest(label)): "
             << "sendVal=" << sendVal << " recvVal_cube=" << recvVal_cube
             << " df/dsendVal=" << computed_deriv << " expected=" << expected_deriv
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // Test 15: AD non-blocking with finishedRequest(Request&) API
    {
        tape.reset();
        tape.setActive();

        int myRank = Pstream::myProcNo();
        int nProcs = Pstream::nProcs();
        int sendTo = (myRank + 1) % nProcs;
        int recvFrom = (myRank - 1 + nProcs) % nProcs;

        scalar sendVal = myRank + 500.0;
        tape.registerInput(sendVal);
        scalar result = sendVal + 10.0;

        Field<scalar> recvBuf(1);
        Field<scalar> sendBuf(1);
        sendBuf[0] = result;

        UPstream::Request recvReq;
        UPstream::Request sendReq;

        UIPstream::read(recvReq, recvFrom, recvBuf.data(), recvBuf.size(), 0, UPstream::worldComm);
        UOPstream::write(sendReq, sendTo, sendBuf.cdata(), sendBuf.size(), 0, UPstream::worldComm);

        // Poll with finishedRequest(Request&)
        while (!UPstream::finishedRequest(recvReq) ||
               !UPstream::finishedRequest(sendReq))
        {
            // Busy wait
        }

        scalar recvResult = recvBuf[0];
        tape.registerOutput(recvResult);
        tape.setPassive();

        // Only master sets gradient to avoid MediPack summing across ranks
        if (Pstream::master())
        {
            recvResult.setGradient(1.0);
        }

        tape.evaluate();

        scalar computed_deriv = sendVal.getGradient();
        scalar expected_deriv = (myRank == nProcs - 1) ? 1.0 : 0.0;

        bool passed = (mag(computed_deriv - expected_deriv) < 1e-10);
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 15 (AD finishedRequest(Request&)): "
             << "sendVal=" << sendVal << " recvResult=" << recvResult
             << " df/dsendVal=" << computed_deriv << " expected=" << expected_deriv
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // Test 16: Non-AD waitRequestPair API
    {
        tape.reset();

        int myRank = Pstream::myProcNo();
        int nProcs = Pstream::nProcs();
        int sendTo = (myRank + 1) % nProcs;
        int recvFrom = (myRank - 1 + nProcs) % nProcs;

        double sendVal1 = myRank + 600.0;
        double sendVal2 = myRank + 700.0;
        double recvVal1 = 0.0;
        double recvVal2 = 0.0;

        // Track request indices manually
        label startIdx = UPstream::nRequests();

        UIPstream::read(UPstream::commsTypes::nonBlocking, recvFrom,
                        reinterpret_cast<char*>(&recvVal1), sizeof(double), 0, UPstream::worldComm);
        UIPstream::read(UPstream::commsTypes::nonBlocking, recvFrom,
                        reinterpret_cast<char*>(&recvVal2), sizeof(double), 1, UPstream::worldComm);
        UOPstream::write(UPstream::commsTypes::nonBlocking, sendTo,
                         reinterpret_cast<const char*>(&sendVal1), sizeof(double), 0, UPstream::worldComm);
        UOPstream::write(UPstream::commsTypes::nonBlocking, sendTo,
                         reinterpret_cast<const char*>(&sendVal2), sizeof(double), 1, UPstream::worldComm);

        // Calculate indices
        label recvReq1 = startIdx;
        label recvReq2 = startIdx + 1;
        label sendReq1 = startIdx + 2;
        label sendReq2 = startIdx + 3;

        // Use waitRequestPair
        UPstream::waitRequestPair(recvReq1, recvReq2);
        UPstream::waitRequestPair(sendReq1, sendReq2);

        bool passed = (recvVal1 == (recvFrom + 600.0)) && (recvVal2 == (recvFrom + 700.0));
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 16 (Non-AD waitRequestPair): "
             << "recvVal1=" << recvVal1 << " recvVal2=" << recvVal2
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // Test 17: Non-AD finishedRequestPair API
    {
        tape.reset();

        int myRank = Pstream::myProcNo();
        int nProcs = Pstream::nProcs();
        int sendTo = (myRank + 1) % nProcs;
        int recvFrom = (myRank - 1 + nProcs) % nProcs;

        double sendVal1 = myRank + 800.0;
        double sendVal2 = myRank + 900.0;
        double recvVal1 = 0.0;
        double recvVal2 = 0.0;

        // Track request indices manually
        label startIdx = UPstream::nRequests();

        UIPstream::read(UPstream::commsTypes::nonBlocking, recvFrom,
                        reinterpret_cast<char*>(&recvVal1), sizeof(double), 0, UPstream::worldComm);
        UIPstream::read(UPstream::commsTypes::nonBlocking, recvFrom,
                        reinterpret_cast<char*>(&recvVal2), sizeof(double), 1, UPstream::worldComm);
        UOPstream::write(UPstream::commsTypes::nonBlocking, sendTo,
                         reinterpret_cast<const char*>(&sendVal1), sizeof(double), 0, UPstream::worldComm);
        UOPstream::write(UPstream::commsTypes::nonBlocking, sendTo,
                         reinterpret_cast<const char*>(&sendVal2), sizeof(double), 1, UPstream::worldComm);

        // Calculate indices
        label recvReq1 = startIdx;
        label recvReq2 = startIdx + 1;
        label sendReq1 = startIdx + 2;
        label sendReq2 = startIdx + 3;

        // Poll with finishedRequestPair
        while (!UPstream::finishedRequestPair(recvReq1, recvReq2) ||
               !UPstream::finishedRequestPair(sendReq1, sendReq2))
        {
            // Busy wait
        }

        bool passed = (recvVal1 == (recvFrom + 800.0)) && (recvVal2 == (recvFrom + 900.0));
        allTestsPassed = allTestsPassed && passed;

        Pout << "Test 17 (Non-AD finishedRequestPair): "
             << "recvVal1=" << recvVal1 << " recvVal2=" << recvVal2
             << " [" << (passed ? "PASS" : "FAIL") << "]" << endl;
    }

    // Print summary and return appropriate exit code
    if (Pstream::master())
    {
        Info << endl << "========================================" << endl;
        if (allTestsPassed)
        {
            Info << "ALL TESTS PASSED" << endl;
        }
        else
        {
            Info << "SOME TESTS FAILED" << endl;
        }
        Info << "========================================" << endl;
    }

    return allTestsPassed ? 0 : 1;
}


// ************************************************************************* //
