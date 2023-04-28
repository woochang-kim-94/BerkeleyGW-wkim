#!/bin/bash
tar --strip=2 -xzf benzothiazole.tgz benzothiazole/1-scf/{QKB,VSC,VKB,WFN_in}
tar -czf data.tgz QKB VKB VSC WFN_in
rm QKB VKB VSC WFN_in
