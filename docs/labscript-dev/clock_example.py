self.clock = [{'start': 0,        'reps': 1,   'step': 1e-3, 'slow_clock_tick':True},
               'WAIT',
              {'start': 1e-3,     'reps': 1,   'step': 1e-6, 'slow_clock_tick':True},
              {'start': 1.001e-3, 'reps': 999, 'step': 1e-6, 'slow_clock_tick':False},
              {'start': 2e-3,     'reps': 1,   'step': 1e-3, 'slow_clock_tick':True}]
