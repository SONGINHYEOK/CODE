
class Robot:
    
    """
    
    [Robot Class]
    Author : 송인혁
    Role : ???
    
    """
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
    
    @staticmethod
    def this_is_robot_class():
        print("yes!!")
    
    #오버라이팅    
    def __str__(self):
        return f"{self.name} robot!!" #<__main__.Robot object at 0x7fd7081dfd60> - > R2-D2 robot!!
     
    def __call__(self):
        print("call!")
        return f"{self.name} call!!"        

droid1 = Robot("R2-D2", 12312321213)
droid1.say_hi()

print(dir(droid1))

print(droid1.__str__())

#python은 모든 것은 객체이다.

droid1()