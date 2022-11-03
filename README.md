# Mhorseshoe

Mhorsehoe algorithm

## 업데이트 및 해야 할 것

approximate algorithm step check 및 최적화 진행

## approximate algorithm 분석

approximate에서 가끔 오류가 발생하는 이유를 찾았다.

총 2가지가 존재한다.

1. eta = 0 인 경우 : 1/eta가 Inf가 된다. 따라서 U matrix에서 Inf, -Inf가 생성된다.
2. epsilon = Inf 인 경우 : eta에서 NA가 샘플링 된다.

eta = 0이 나오는 경우는, beta값이 너무 큰 경우에 해당한다.

eta = NA가 나오는 경우는 beta값이 엄청커서 Inf에 다가설 정도에 해당한다.

eta * xi의 최소값으로 10^(-30)을 설정했더니, beta가 엄청나게 커졌다. 아마 diagonal_delta 값이 너무 커도 beta값이 뭉개지는 모양이다.

그러므로, eta * xi == 0 인 경우 truncate를 해서 최소값을 정해줘야 한다. 나는 10^(-8)로 정했다. 주변 non-zero의 eta * xi 값과 큰 차이가 나지 않으면서도 active set임을 증명할 정도이면 될 것 같다.


정말 본질적인 원인은, eta가 0이 되는게 아니다.

오히려 eta가 딱 0이 되어주면, active set에서 xi의 값에 관계없이 당연히 통과할 수 있는 보증 수표가 된다.

가장 큰 원인은 eta = 0인 부분이 있을 때, diagonal이 0이 되는 부분이 생기는데,

이 부분을 해결하지 않을 경우, 일단 diagonal_delta값이 Inf가 되는 부분이 생기고, 이러면 U matrix에서 반드시 Inf, -Inf가 특정 행에서 생성된다.

그럼 Inf-Inf = NaN이므로 U %*% v_star가 NaN으로 도출되므로 new_beta가 NaN으로 도출된다.

그래서 diagonal에서 



## git 동기화 관련

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
5. git push origin --delete brach.name을 할 경우 repository의 branch를 삭제할 수 있다.
