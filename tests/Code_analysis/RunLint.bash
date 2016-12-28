#!/bin/bash

./cpplint.py --filter=-build/include,-whitespace,-legal ../../src/* > LintTest.log 2>&1
sed -i '/Ignoring /d' LintTest.log
