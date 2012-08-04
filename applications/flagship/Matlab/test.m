clear

matDiagScalarSegregatedCpu
matDiagScalarSegregatedGpu
matDiagBlockSegregatedCpu
matDiagBlockSegregatedGpu

disp("DIAG: CPU: scalar - GPU scalar")
matDiagScalarSegregatedCpu-matDiagScalarSegregatedGpu

disp("DIAG: CPU: block - GPU block")
matDiagBlockSegregatedCpu-matDiagBlockSegregatedGpu


matDiagScalarCoupledCpu
matDiagScalarCoupledGpu
matDiagBlockCoupledCpu
matDiagBlockCoupledGpu

disp("DIAG: CPU: scalar - GPU scalar")
matDiagScalarCoupledCpu-matDiagScalarCoupledGpu

disp("DIAG: CPU: block - GPU block")
matDiagBlockCoupledCpu-matDiagBlockCoupledGpu

matFullScalarSegregatedCpu
matFullScalarSegregatedGpu
matFullBlockSegregatedCpu
matFullBlockSegregatedGpu

disp("FULL: CPU: scalar - GPU scalar")
matFullScalarSegregatedCpu-matFullScalarSegregatedGpu

disp("FULL: CPU: block - GPU block")
matFullBlockSegregatedCpu-matFullBlockSegregatedGpu


matFullScalarCoupledCpu
matFullScalarCoupledGpu
matFullBlockCoupledCpu
matFullBlockCoupledGpu

disp("FULL: CPU: scalar - GPU scalar")
matFullScalarCoupledCpu-matFullScalarCoupledGpu

disp("FULL: CPU: block - GPU block")
matFullBlockCoupledCpu-matFullBlockCoupledGpu

