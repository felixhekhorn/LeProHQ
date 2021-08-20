Using [cProfile](https://docs.python.org/3/library/profile.html) and [flameprof](https://github.com/baverman/flameprof) you can check the performance:

```
python -m cProfile -o myscript.prof integration.py 
flameprof myscript.prof > output.svg
```
