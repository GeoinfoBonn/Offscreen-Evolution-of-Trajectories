# Offscreen-Evolution-of-Trajectories

This repository contains the source code accompanying the publication *Visualizing the Off-Screen Evolution of Trajectories* in KN by Axel Forsch, Friederike Amann and Jan-Henrik Haunert.

## Dependencies

* Java 11 JRE
* Maven

## Building the project

The project uses Maven as dependency management tool. The project can be compiled using the following comand in the repository root directory:

```sh
mvn clean package
```

The compiled jar-file can be found in the build directory `target/` with the name `OffscreenEvolution.jar`.

## Running the tool

The tool can be run from the command line using flags to specify the input parameters. The following is an example call using the input data provided in `src/main/resources/`:

```sh
java -jar target/OffscreenEvolution.jar -in src/main/resources/tracks/allTracks.shp -bounds 797079.042,6582411.495,806112.618,6588239.608 -out out/ -beta 0.1 -kappa .90 -t 150 -ns 3 -ol 0.5 -vmax 7 -var ALL_POINTS
```

A list of all flags with a short description can be found by calling:

```sh
java -jar target/OffscreenEvolution.jar -h`
```
