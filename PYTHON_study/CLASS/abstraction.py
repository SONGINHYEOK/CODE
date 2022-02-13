"""

#* 추상화 : abstraction
#* 불필요한 정보는 숨기고 중요한(필요한) 정보만능ㄹ 표현함으로써
#* 공통의 속성 값이나 행위(methods)를 하나로 묶어 이름을 붙이는 것이다.

"""
#class로 구성
class Robot:
    
    #class 변수 : 인스턴스들이 공유하는 변수
    population = 0
    
    #생성자 함수
    def __init__(self, name, code):
        self.name = name #인스턴스 변수
        self.code = code #인스턴스 변수
        Robot.population +=1
    
    #인스턴스 메서드
    def say_hi(self):
        #code...
        print(f"Greetings, my masters call me {self.name}.")
    
    #인스턴스 메서드    
    def cal_add(self, a, b):
        return a + b
    
    #인스턴스 메서드
    def die(self):
        print(f"{self.name} is being destroyed!")
        Robot.population -= 1
        if Robot.population == 0:
            print(f"{self.name} was the last one.")
        else:
            print(f"There are still {Robot.population} robots workig.")
    
    #class method    
    @classmethod
    def how_many(cls):
        print(f"We have {cls.population} robots.")    

print(Robot.population)

siri = Robot("siri", 123124511515)
jarvis = Robot("jarvis", 12355512134)
bixby = Robot("bixby", 125901923213)

print(siri.name)

siri.say_hi()
siri.cal_add(1,3)

Robot.how_many()




#개별로 구성할 때
"""
siri_name = "siri"

siri_code = 123124511515

def siri_say_hi():
    #code ....
    print("say hello!!, my name is siri")
    
def siri_add_cal():
    return 2 + 3

def siri_die():
    #code ...
    print("siri die")
    
    
jarvis_name = 'javis' 
jarvis_code = 12355512134

#.....

""" 