from functools import partial
import time
from concurrent.futures import ThreadPoolExecutor
from tornado import gen
from bokeh.document import without_document_lock
from bokeh.models import ColumnDataSource
from bokeh.plotting import curdoc, figure
import numpy as np

def add_noise():
  return (np.random.choice([-1,1])*np.random.rand())/10

x=np.arange(0,1,0.1)

def cos_n(x):
    y=[y+add_noise() for y in np.cos(x)]
    return np.array(y)    

source = ColumnDataSource(data=dict(x=x, y=np.sin(x), z=cos_n(x)))

i = 0

doc = curdoc()

max_x=20

executor = ThreadPoolExecutor(max_workers=2)
plot = figure(x_range=[source.data['x'][-1], source.data['x'][-1]+max_x], y_range=[-1.3,1.3])
plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6, color='blue', legend_label="sin(x)")
plot.line('x', 'z', source=source, line_width=3, line_alpha=0.6, color='red', legend_label="cos(x)")

# Funcion con la politica de bloque de tareas
def blocking_task(i):
    time.sleep(1)
    return i

#Esta funcion es necesaria para la funcion que actualiza unlocked actualice
#la información de forma segura (debe nombrarse locked_update)
@gen.coroutine
def locked_update(i):
    last=source.data['x'][-1]
    x=[xi+add_noise() for xi in np.arange(last+0.1,last+0.5,0.1)] 
    source.stream(dict(x=x,y=np.sin(x), z=cos_n(x)))
    plot.x_range.end = max(source.data['x'][-1],max_x)
    plot.x_range.start = plot.x_range.end-max_x

# La funcion unlocked que no bloquea ejecuaciones de otras llamada
# mientras se esta ejecutando 
@gen.coroutine
@without_document_lock
def unlocked_task():
    global i
    i += 1
    res = yield executor.submit(blocking_task, i)
    doc.add_next_tick_callback(partial(locked_update, i=res))

# Este es el callback sin bloqueo
@gen.coroutine
def update():
    last=source.data['x'][-1]
    x=[xi+add_noise() for xi in np.arange(last+0.1,last+0.5,0.1)]	
    source.stream(dict(x=x,y=np.sin(x), z=cos_n(x)))
    plot.x_range.end = max(source.data['x'][-1],max_x)
    plot.x_range.start = plot.x_range.end-max_x


doc.add_periodic_callback(unlocked_task, 1000)
doc.add_periodic_callback(update, 300)
doc.add_root(plot)

