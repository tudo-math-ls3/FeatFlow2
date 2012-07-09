@echo off

echo.
echo Cleaning up intermediate build files...

rem Delete Visual C++ user-config files; these are created automatically
del /q      visual_studio\*.vcxproj.user

rem Delete Visual Fortran user-config files; these are created automatically
del /q /a:h visual_studio\*.u2d

rem Delete object and other intermediate files generated during build
rmdir /s /q visual_studio\obj

echo Cleanup finished