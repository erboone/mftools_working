import unittest

class test: 
    foo = 'this works'

    def __init__(self):
        self.bar = 'this too'

class newtest(test):
    pass


if __name__ == "__main__":
    var = 'string'

    t = newtest()
    
    print(t.foo)
    print(t.bar)


