#!/bin/bash
sed '/>/a %' $1|sed -z 's/\n//g'|sed 's/>/\n>/g'|grep -v '%.*[XNYRKMWSBHVD].*'|sed 's/%/\n/g'
