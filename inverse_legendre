import time, math

def prime_factors_sieve(limit):
    result = [{} for _ in range(limit + 1)]
    for i in range(2, limit + 1):
        if len(result[i]) == 0:
            #prime number found
            for j in range(i, limit + 1, i):
                n = j
                if i in result[j]:
                    while n % i == 0:
                        n //= i
                        result[j][i] += 1
                else:
                    result[j][i] = 1
                    n //= i
                    while n % i == 0:
                        n //= i
                        result[j][i] += 1
    return result

def legendre_factorial(n, p):
    #Modified legendre factorial function, finds number of factors of p in n!
    total = 0
    for i in range(1, int(math.floor(math.log(n, p))) + 1):
        total += n // (p ** i)
    return total

def inv_legendre_factorial(p, e):
    #Finds smallest number, n, such that p^e divides n!
    #Use bisection method
    low = p #Low guess is n = p
    high = p*e #high guess is n = p*e
    while high - low > 1:
        mid = (low + high)//2
        if legendre_factorial(mid, p) >= e:
            #If middle number has more factors of p than needed
            #make high = mid
            high = mid
        else:
            #Otherwise low = mid
            low = mid
    #After high and low are one integer apart we test high and low
    if legendre_factorial(low, p) >= e:
        #If low has enough prime factors, we return low
        return low
    else:
        #Otherwise return high
        return high
