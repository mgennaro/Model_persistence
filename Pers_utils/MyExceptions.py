class CustExc1M1V(Exception):
    '''
      Custom exception class to raise exceptions with one messgae and one value
    '''

    def __init__(self, msg, val):
        self.msg = msg
        self.val = val


class CustExc1M3V(Exception):
    '''
      Custom exception class to raise exceptions with one messgae and 3 values
    '''

    def __init__(self, msg, val1, val2, val3):
        self.msg = msg
        self.val1 = val1
        self.val2 = val2
        self.val3 = val3
