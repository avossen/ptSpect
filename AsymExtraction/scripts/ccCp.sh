#/bin/bash
a=$1
destination=${a##*/}
echo destination: $destination
echo "ssh -L 3128:192.168.1.130:3128 vossen@sshcc2.kek.jp  -t scp vossen@login.cc.kek.jp:~/$1  /tmp/ && scp vossen@sshcc2.kek.jp:/tmp/$destination"
ssh vossen@sshcc2.kek.jp  -t rm "/tmp/$destination"
ssh -A -L 3128:192.168.1.130:3128 vossen@sshcc2.kek.jp  -t scp "vossen@login.cc.kek.jp:~/$1  /tmp/"
scp "vossen@sshcc2.kek.jp:/tmp/$destination" .
ssh vossen@sshcc2.kek.jp  -t rm "/tmp/$destination"
