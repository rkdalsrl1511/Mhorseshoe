# Mhorseshoe algorithm

## 간단 사용법


```r
# 시뮬레이션 데이터 생성
N <- 200
p <- 1000
W <- Mhorseshoe::make_W(N = N, p = p)

# design matrix 표준화 : 반드시 이 함수로 해야함.
standardized_W <- scale(W[,-1], center = TRUE, scale = TRUE)

# 시뮬레이션 데이터 response variable 만들기
# non_zero : true signal의 수 설정
# SD : linear model에서 error term의 표준편차
# fixed_coefficients : 1로 설정할 경우 모든 true signal의 coefficients가 1로 설정 됨
# 설정하지 않는 경우, -100~100 사이의 값으로 랜덤하게 할당 됨
# c(1,2,3,4,5,6,7,8,9,10)으로 입력할 경우 각각 1,2,3,4,5,6,7,8,9,10으로 설정 됨
response_list <- Mhorseshoe::make_response(standardized_W, non_zero = 10, 
                                           SD = 1, 
                                           fixed_coefficients = 1)

# response variable
z <- response_list[[1]]

# true signal의 열 위치
non_zero_index <- response_list[[2]]

# true signal의 coefficients
non_zero_coefficents <- response_list[[3]]

# approximate 알고리즘 데이터에 적합(모든 설정을 디폴트로)
fit_amcmc <- Mhorseshoe::approximate_algorithm(standardized_W, z,
                                               iteration = 10000)

# 수정한 알고리즘 데이터에 적합(모든 설정을 디폴트로)
fit_Mamcmc <- Mhorseshoe::modified_approximate_algorithm(standardized_W, z,
                                                         iteration = 10000)

```

## 업데이트

### 11/09

#### modified algorithm 수정

meff 관련, xi값 수정






## 해야 할 것 

1. truncate rejection sampler

2. meff 관련 수정



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


깃허브의 private repository이므로

https://rkdalsrl1511@github.com/rkdalsrl1511/Mhorseshoe.git

으로 입력해야 오류가 발생하지 않음

# approximate algorithm

approximate algorithm의 딜레마

eta를 truncate할 경우, 작은 coefficients가 0으로 shrinkage 해버린다.

eta를 truncate하지 않는 경우, eta = 0이 되어서 알고리즘이 망가질 수 있다.

또한, truncate를 얼마나 하느냐도 매우 중요하다.

너무 작은 값으로 하면, 그럼 또 error가 발생하고, 

너무 큰 값으로 하면 xi는 성장하지 않아서, shrinkage가 과도하게 발생한다.
