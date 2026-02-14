# meltPoolFoam

OpenFOAM solver for simulation melt-pool dynamics.
=======

OpenFOAM solver for simulation melt-pool dynamics.

## Development Environment Setup (VS Code)

This project uses **Visual Studio Code** with a **Dev Container** to provide a consistent, pre-configured OpenFOAM development environment. This ensures that all dependencies, compilers, and code intelligence (powered by `ccls`) work automatically without manual configuration.

### Prerequisites

1.  **Docker Desktop** (running and up to date)
2.  **Visual Studio Code**
3.  **Dev Containers Extension** for VS Code

### macOS Users (Critical: Case Sensitivity)

By default, the macOS file system (APFS) is **case-insensitive**, which causes name conflicts in OpenFOAM (e.g., it cannot distinguish between `scalar` and `Scalar` files).

To fix this and ensure correct behavior, you **must** create a case-sensitive volume (disk image) and store your OpenFOAM project there.

Run the following commands in your terminal **before** opening the project:

```bash
# 1. Install the specialized macOS file system driver for OpenFOAM
sudo curl -o /usr/local/bin/openfoam-macos-file-system http://dl.openfoam.org/docker/openfoam-macos-file-system
sudo chmod 755 /usr/local/bin/openfoam-macos-file-system

# 2. Create a case-sensitive volume named 'my_openfoam'
sudo openfoam-macos-file-system -v my_openfoam create

# 3. Mount the volume locally
mkdir -p my_openfoam
sudo openfoam-macos-file-system -v my_openfoam mount
