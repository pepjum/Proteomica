#!/bin/bash

find /mnt/beegfs/spack/opt/spack/linux-centos7-x86_64/gcc-4.8.5/ -mindepth 2 -type d -iname "bin" | while read line; do
    cp -n "$line"/* /mnt/beegfs/agarin/bin 2>/dev/null
done

find /mnt/beegfs/spack/opt/spack/linux-centos7-x86_64/gcc-4.9.4/ -mindepth 2 -type d -iname "bin" | while read line; do
    cp "$line"/* /mnt/beegfs/agarin/bin 2>/dev/null
done