#!/bin/bash
tar --strip=2 -xzf Zn2O2.tgz Zn2O2/{1-scf/QKB,1-scf/VSC,2.1-nscf/VKB,2.1-nscf/WFN_in}
tar -czf data.tgz QKB VKB VSC WFN_in
rm QKB VKB VSC WFN_in
