#!/bin/bash

git archive --prefix=forqs/ master | bzip2 > forqs.tbz
scp forqs.tbz dkessner@hoffman2.idre.ucla.edu:~/dev

