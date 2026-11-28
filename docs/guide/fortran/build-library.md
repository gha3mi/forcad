title: Build the Fortran library
category: Setup

Install [Pixi](https://pixi.sh/) and Git. On macOS, install the Xcode command-line tools. If you already have the checkout, start in its root folder.

```sh
git clone --branch dev https://github.com/gha3mi/forcad.git
cd forcad
pixi run -e gnu build --profile release
```

The `gnu` environment supplies GNU Fortran. fpm reports a successful build or that the project is up to date. The [build reference](../reference/build.html) lists other compiler environments and the debug test command.
