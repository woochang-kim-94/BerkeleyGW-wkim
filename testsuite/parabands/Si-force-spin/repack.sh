#!/bin/bash
#tar --strip=2 -xzf Si-force-spin.tgz Si-force-spin/{1-scf/QKB,1-scf/VSC,2.1-nscf-no-cons/VKB,2.1-nscf-no-cons/WFN_in}
tar --strip=2 -xzf Si-force-spin.tgz Si-force-spin/2.1-nscf-no-cons/{QKB,VSC,VKB,WFN_in}
tar -czf data.tgz QKB VKB VSC WFN_in
rm QKB VKB VSC WFN_in
