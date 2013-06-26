function rrr_manipulator

    tSpan = linspace(0, 5, 100); % seconds
    initialConditions = [pi/4; 0; pi/10; 0; 0; 0];

    p.m1 = 1; % kg
    p.m2 = 1; % kg
    p.m3 = 1; % kg
    p.L1 = 1; % kg
    p.L2 = 1; % kg
    p.L3 = 1; % kg
    p.g = 9.81; % m/s

    [times states] = ode45(@(t, x) equationsOfMotion(t, x, p), ...
        tSpan, initialConditions);

    Q1s = states(:, 1);
    Q2s = states(:, 2);
    Q3s = states(:, 3);

    X1s = p.L1 * cos(Q1s);
    Y1s = p.L1 * sin(Q1s);

    X2s = X1s + p.L2 * cos(Q2s);
    Y2s = Y1s + p.L2 * sin(Q2s);

    X3s = X2s + p.L3 * cos(Q3s);
    Y3s = Y2s + p.L3 * sin(Q3s);

    hf = figure;
    hold on;
    axis equal;
    for iTime = 1:length(times)
        cla;
        axis([-3 3 -3 3]);
        plot([0, X1s(iTime), X2s(iTime), X3s(iTime)], ...
             [0, Y1s(iTime), Y2s(iTime), Y3s(iTime)], '-o');
        pause(0.05);
    end

end

function tau = control(t, x, p)



    % TODO
    % first: gravity compensation.
    %tau = Jacobian' * F;

end

function xDot = equationsOfMotion(t, x, p)

    q1 = x(1);
    q2 = x(2);
    q3 = x(3);
    u1 = x(4);
    u2 = x(5);
    u3 = x(6);

    q1d = u1;
    q2d = u2;
    q3d = u3;

    m1 = p.m1;
    m2 = p.m2;
    m3 = p.m3;
    L1 = p.L1;
    L2 = p.L2;
    L3 = p.L3;

    m23 = m2 + m3;
    m123 = m1 + m23;

    MM11 = m123 * L1^2;
    MM22 = m23 * L2^2;
    MM33 = m3 * L3^2;

    MM12 = m23 * L1 * L2 * cos(q2 - q1);
    MM13 = m3 * L1 * L3 * cos(q3 - q1);
    MM23 = m3 * L2 * L3 * cos(q3 - q2);

    % Dynamics.
    q12 = q1 + q2;
    q123 = q1 + q2 + q3;

    J11 = -p.L1*sin(q1) - p.L2*sin(q12) - p.L3*sin(q123);
    J12 = -p.L2*sin(q12) - p.L3*sin(q123);
    J13 = -p.L3*sin(q123);
    J21 = p.L1*cos(q1) + p.L2*cos(q12) + p.L3*cos(q123);
    J22 = p.L2*cos(q12) + p.L3*cos(q123);
    J23 = p.L3*cos(q123);
    Jacobian = [J11 J12 J13;
                J21 J22 J23];

    MM11 = m1*L1^2 + m2*(L1^2+L2^2+2*L1*L2*cos(q2)) + ...
        m3*(L1^2+L2^2+L3^2+2*L1*L2*cos(q2)+2*L2*L3*cos(q3)+2*L1*L3*cos(q2+q3));
    MM12 = L2*m2*(L2+L1*cos(q2)) + ...
        m3*(L2^2+L3^2+L1*L2*cos(q2)+2*L2*L3*cos(q3)+L1*L3*cos(q2+q3));
    MM13 = L3*m3*(L3+L2*cos(q3)+L1*cos(q2+q3));
    MM22 = m2*L2^2 + m3*(L2^2+L3^2+2*L2*L3*cos(q3));
    MM23 = L3*m3*(L3+L2*cos(q3));
    MM33 = m3*L3^2;

    MassMatrix = [MM11 MM12 MM13;
                  MM12 MM22 MM23;
                  MM13 MM23 MM33];

    QV1 = L1*L2*m2*sin(q2)*(q1d^2-(q1d+q2d)^2) + ...
        m3*(L1*L2*sin(q2)*q1d^2+L1*L3*sin(q2+q3)*q1d^2+...
        L2*L3*sin(q3)*(q1d+q2d)^2-L1*L2*sin(q2)*(q1d+q2d)^2-...
        L2*L3*sin(q3)*(q1d+q2d+q3d)^2-L1*L3*sin(q2+q3)*(q1d+q2d+q3d)^2);

    QV2 = L1*L2*m2*sin(q2)*q1d^2 + ...
        m3*(L1*L2*sin(q2)*q1d^2+L1*L3*sin(q2+q3)*q1d^2+...
        L2*L3*sin(q3)*(q1d+q2d)^2-L2*L3*sin(q3)*(q1d+q2d+q3d)^2);

    QV3 = L3*m3*(L1*sin(q2+q3)*q1d^2+L2*sin(q3)*(q1d+q2d)^2);

    QuadraticVelocity = [QV1;
                         QV2;
                         QV3];

    G1 = L1*m1*cos(q1) + m2*(L1*cos(q1)+L2*cos(q1+q2)) + ...
        m3*(L1*cos(q1)+L2*cos(q1+q2)+L3*cos(q1+q2+q3));
    G2 = L2*m2*cos(q1+q2) + m3*(L2*cos(q1+q2)+L3*cos(q1+q2+q3));
    G3 = L3*m3*cos(q1+q2+q3);
    Gravity = p.g * [G1;
                     G2;
                     G3];


    invMassMatrix = inv(MassMatrix);
    TaskSpaceMassMatrix = inv(Jacobian * invMassMatrix * Jacobian');

    DynConsistentJacInverse = invMassMatrix * Jacobian' * TaskSpaceMassMatrix;
    %JacPseudoinverse = Jacobian' * inv(Jacobian * Jacobian');

    TaskSpaceGravity = DynConsistentJacInverse' * Gravity;

    EndEffectorForce = TaskSpaceGravity;

    % TODO Nullspace action.
    Nullspace = eye(3) - DynConsistentJacInverse * Jacobian;

    Actuators = Jacobian' * EndEffectorForce + Nullspace' * ;


    % Solve.
    qDoubleDot = MassMatrix \ (-QuadraticVelocity - Gravity + Actuators);

    u1d = qDoubleDot(1);
    u2d = qDoubleDot(2);
    u3d = qDoubleDot(3);

    xDot = [q1d; q2d; q3d; u1d; u2d; u3d];

end



