# Mhorseshoe

Mhorsehoe algorithm

## git 동기화

우선 git을 설치하고 셸에서 다음을 입력한다.

git config --global user.name ""
git config --global user.email ""

이후 R 프로젝트에서 project options -> version control system을 None에서 git으로 바꾼다.

자동으로 재시작이 되며, 이제 git에서 repository를 새롭게 만들어준다(프로젝트와 동일한 이름으로).

그리고 origin을 설정하여 git의 repository와 프로젝트를 연동시켜야 한다.

1. 셸에 git remote add origin https://github.com/rkdalsrl1511/Mhorseshoe.git 입력
2. 셸에 git push -u origin master 입력
3. default branch는 settings에서 바꿀 수 있음.
4. git remote remove origin을 할 경우, 설정된 origin을 초기화할 수 있다.
