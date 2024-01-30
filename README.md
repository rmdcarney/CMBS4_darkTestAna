# Dark tests analysis

Transforms, fits, and plots to analyse dark test data.

## Compile & run

Use makefile. `clean` target exists.
Run like: `./bin/programName $inputFile $outputPrefix `

## Structure

Main executable src in `exec` dir, `src` and `include` contain namespace body & header, `scripts` contains bash scripts. 
Compiled executable in `bin`.
3 namespaces: utils, rootUtils, num

