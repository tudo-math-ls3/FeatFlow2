This "libraries/kernel" subdirectory contains a Microsoft Visual Studio
project for creating a library from the kernel/ source files.
It's not used within UNIX MAKE scripts but simplifies the handling
of the kernel sourcefiles within Projects in Windows.

The library build with the project files here is associated
to any Windows application. As the Visual Studio always checks
also library source files whether they are changed or not, changing a
source file in the Widows-kernel-library will automatically lead
to a rebuild of the library.
