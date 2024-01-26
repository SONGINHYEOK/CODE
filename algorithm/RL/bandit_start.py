import numpy as np

class Bandit:
    def __init__(self, arms=10): #arms 슬롯머신 대수
        self.rates = np.random.rand(arms) #슬롯머신 각각의 승률 설정(무작위)
    
    def play(self, arm):
        rate = self.rates[arm] 
        if rate > np.random.rand():
            return 1
        else:
            return 0
        
bandit = Bandit()

Qs = np.zeros(10)
ns = np.zeros(10)


for n in range(10):
    action = np.random.randint(0, 10)
    reward = bandit.play(action)
    
    ns[action] += 1 #action번째 슬롯머신을 플레이한 횟수 증가
    Qs[action] += (reward - Qs[action] / ns[action])
    print(Qs)