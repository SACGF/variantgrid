#!/bin/sh

find /tau/data/clinical/unaligned -maxdepth 3 -name "RTAComplete.txt" -exec dirname {} \;

