# Mhorseshoe algorithm

## 간단 사용법


```r
# 시뮬레이션 데이터 생성
N <- 200
p <- 1000
W <- Mhorseshoe::make_W(N = N, p = p)

############################ make_response function #############################
# non_zero : nonzero parameter의 수 설정                                        #
# SD : linear regression model에서 error term의 표준편차                        #
# fixed_coefficients : NULL인 경우, -100~100사이의 값 랜덤 할당. 숫자를 입력하면#
# 해당 숫자를 nonzero coefficients에 할당. vector form도 가능.                  #
# c(1,2,3,4,5,6,7,8,9,10)으로 입력할 경우 각각 1,2,3,4,5,6,7,8,9,10으로 설정 됨 #
# 결과는 list 형태를 반환                                                       #
#################################################################################

response_list <- Mhorseshoe::make_response(standardized_W,
                                           non_zero = 10, 
                                           SD = 1, 
                                           fixed_coefficients = 5)

# response variable
z <- response_list[[1]]

# true signal의 index
non_zero_index <- response_list[[2]]

# true signal의 coefficients
non_zero_coefficents <- response_list[[3]]

# default 설정으로 approximate 알고리즘 데이터에 적합
fit_amcmc <- Mhorseshoe::approximate_algorithm(standardized_W, z,
                                               iteration = 10000)

# default 설정으로 수정한 알고리즘 데이터에 적합
fit_Mamcmc <- Mhorseshoe::modified_approximate_algorithm(standardized_W, z,
                                                         iteration = 10000)
```

## 업데이트

### 11/09

- modified algorithm 수정
- meff 계산법 수정, xi값 수정

### 11/10

- modified algorithm 코드 수정 : threshold값 설정할 필요 없도록 수정
- truncate rejection sampler 만들었다가 효과가 없음을 확인

### 12/05

- rejection sampler 수정(Epsilon 값에 따른 eta, xi값 변동 확인)

## 해야 할 것


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
