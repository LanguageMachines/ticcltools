#!/bin/sh

if [ -n "$*" ]; then
    exec "$@"
else
    sh
fi
