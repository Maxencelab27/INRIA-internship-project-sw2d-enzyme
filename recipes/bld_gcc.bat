:: rmdir /s /q build
mkdir build
cd build

:: set build configuration to command-line argument value if it exists, else to Release
IF "%1"=="" (set BUILD_CONFIG=RELEASE) ELSE (set BUILD_CONFIG=%1)


cmake .. -G "MinGW Makefiles" ^
    -DCMAKE_BUILD_TYPE=%BUILD_CONFIG% ^
    -DCMAKE_MAKE_PROGRAM=mingw32-make ^
    -DCMAKE_CXX_COMPILER=g++

if errorlevel 1 exit 1

cmake --build . --config %BUILD_CONFIG%
if errorlevel 1 exit 1

cd ..