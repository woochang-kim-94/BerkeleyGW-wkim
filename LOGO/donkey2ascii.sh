#!/bin/bash

jp2a --width=81 --grayscale -i --chars=' ...ooxx@@' donkey.jpg | sed -e "s/\[37m//g;s/\[0m//g;s/ $/'/;s/^ /write(6,'(a)') '/" > donkey.f90
