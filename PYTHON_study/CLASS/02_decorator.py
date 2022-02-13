# decirator

def copyright(func):
    
    def new_func():
        print("@ INHYEOK SONG")
        func()
        
    return new_func

@copyright
def smile():
    print("^_^")
    
    
#smile = copyright(smile)

smile()