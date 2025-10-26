
rm -f cadenas_random.mod
rm -f -r polymer.x*

#flags="-O -Wall -g -fcheck=all -fbacktrace -ffpe-trap=overflow"
#flags="-O -Wall -g -fcheck=all -fbacktrace"

flags="-O3 -Wall"

gfortran $flags cadenas_random.f95 aux-main.f95 -o polymer.x

rm -f cadenas_random.mod

