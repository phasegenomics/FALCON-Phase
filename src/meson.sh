mkdir -p build-meson
meson --buildtype=release ./build-meson
#meson --buildtype=debugoptimized ./build-meson
ninja -C build-meson -v

# executable is now in ./build-meson/falcon-phase
