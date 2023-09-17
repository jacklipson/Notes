import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from numpy import exp, cos, arccos

# used starter code for sliders and button from matplotlib site + stackoverflow for debug but all work and equations are mine!

# defining x(t)
def f(t, beta, omega):
    c = (beta**2 - omega**2)
    if c > 0: # overdamping
        d = c**.5
        return ((beta + d) * exp(t * (-beta + d)) + (-beta + d) * exp(t * (-beta - d))) / (2 * d)
    elif c == 0: #critical damping
        return exp(-beta * t) * (1 + beta * t)
    else: # under damping
        d = (-c)**.5
        return exp(-beta * t) * cos(d * t)

# refining plot w/ 1000 intervals from 0 to 25s
t = np.linspace(0, 25, 1000)

# default parameters
init_beta = 1.0
init_omega = 1.0

# figure plot and bounds
fig, ax = plt.subplots()
line, = ax.plot(t, f(t, init_beta, init_omega), lw=2)
ax.set_xlabel("time (s)")
ax.set_ylabel("x(t) (m)")
plt.title("A damped oscillator")

ax.set_xlim([0,25])
ax.set_ylim([-1,1])

# make room for sliders
fig.subplots_adjust(left=0.25, bottom=0.25)

# make sliders.
axbeta = fig.add_axes([0.25, 0.1, 0.65, 0.03])
beta_slider = Slider(
    ax=axbeta,
    label=r'$\beta = \frac{b}{2m}}$',
    valmin=0,
    valmax=10,
    valinit=init_beta,
)
axomega = fig.add_axes([0.25, 0.05, 0.65, 0.03])
omega_slider = Slider(
    ax=axomega,
    label=r'$\omega_0 = \sqrt{\frac{k}{m}}$',
    valmin=0,
    valmax=10,
    valinit=init_omega,
)

# reset plot data when parameter changes
def update(val):
    line.set_ydata(f(t, beta_slider.val, omega_slider.val))
    fig.canvas.draw_idle()

# update sliders
beta_slider.on_changed(update)
omega_slider.on_changed(update)

# make button to reset sliders
resetax = fig.add_axes([0.02, 0.095, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')

def reset(event):
    beta_slider.reset()
    omega_slider.reset()
button.on_clicked(reset)

plt.show()