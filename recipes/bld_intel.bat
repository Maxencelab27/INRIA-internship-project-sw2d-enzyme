:: rmdir /s /q build
mkdir build
cd build

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"

:: set build configuration to command-line argument value if it exists, else to Release
IF "%1"=="" (set BUILD_CONFIG=RELEASE) ELSE (set BUILD_CONFIG=%1)


cmake .. -G Ninja ^
    -DCMAKE_BUILD_TYPE=%BUILD_CONFIG% ^
    -DCMAKE_CXX_COMPILER=icx ^
    -DCMAKE_MAKE_PROGRAM=ninja

if errorlevel 1 exit 1

cmake --build . --config %BUILD_CONFIG%
if errorlevel 1 exit 1

cd ..
