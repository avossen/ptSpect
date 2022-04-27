#/bin/bash
a=$1
destination=${a##*/}
echo destination: $destination
echo "ssh -L 3128:192.168.1.130:3128 vossen@sshcc1.kek.jp  -t scp vossen@login.cc.kek.jp:~/$1 -r /tmp/ && scp -r vossen@sshcc1.kek.jp:/tmp/$destination"
ssh vossen@sshcc1.kek.jp  -t rm -rf "/tmp/$destination/"
ssh -A -L 3128:192.168.1.130:3128 vossen@sshcc1.kek.jp  -t scp -r "vossen@login.cc.kek.jp:~/$1  /tmp/"
scp -r "vossen@sshcc1.kek.jp:/tmp/$destination" "./$destination/"
ssh vossen@sshcc1.kek.jp  -t rm -rf "/tmp/$destination/"
