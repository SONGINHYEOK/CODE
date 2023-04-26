with open ("/Users/song-inhyeok/CODING/PROTOTYPE/BIChem/NT006/92_final_exescf.txt", 'r') as f:
    lines=f.readlines()

    result1 = set(lines)
    #print(f"set(arr)      : {result1}")

    result2 = list(result1)  # list(set(arr))
    #print(f"list(set(arr) : {result2}")
   
    lines = [line.rstrip('\n') for line in result2]
    print(lines)