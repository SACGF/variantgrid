#!/bin/sh

# Only look in idt_* dirs - as that's where new stuff is put
find -L /tau/data/clinical/unaligned/idt_* -maxdepth 2 -name "RTAComplete.txt" -exec dirname {} \;

