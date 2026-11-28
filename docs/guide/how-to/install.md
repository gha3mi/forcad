title: Launch ForCAD Studio
category: Getting started

<span id="launch-studio"></span>

Install [Pixi](https://pixi.sh/) and Git. On macOS, install the Xcode command-line tools. On Windows, install the Visual Studio Build Tools C++ workload and open a terminal where the tools are available.

```sh
git clone --branch dev https://github.com/gha3mi/forcad.git
cd forcad
pixi run start
```

If you already have the checkout, run `pixi run start` from its root. The first run installs dependencies and builds the application. When Home appears, select **Chess** to follow [Explore a project in Studio](../tutorials/explore-studio.html).

<span id="build-and-test-the-fortran-library"></span><span id="run-the-bundled-curve-example"></span>
For library development, see [Build the Fortran library](../fortran/build-library.html) or [Run the bundled curve example](../fortran/run-curve-example.html).
