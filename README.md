# numsim

## Build
Either install vtk with something like `sudo apt install libvtk7-dev libvtk7.1`, or (since the equivalent `sudo pacman -S vtk` seems to miss some dependencies on arch and therefore makes a build impossible) use vcpkg:
```
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
./bootstrap-vcpkg.sh (you can disable telemetry here)
./vcpkg install vtk
```
Clone numsim with `git@github.com:julianschmoll/numsim.git`, then proceed to build:

```
cd numsim
rm -rf build
cd build
cmake .. (-DCMAKE_TOOLCHAIN_FILE=~/vcpkg/scripts/buildsystems/vcpkg.cmake)
make
```
