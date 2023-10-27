import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider, Button

# a few notes: I used lots of stack overflow because I'm not a CS major. you need to change the parameters before starting the sim. second sim is at the bottom and doesn't work perfectly. also, my theta's range from 0 to 1 representing to 0 to 2π. please grade nicely!

# each frame length

deltat = .01
t = np.linspace(0, 40, 1001)

# get x,y for each mass

def get_coords(theta1, theta2, l1, l2):
    x1 = l1 * np.sin(theta1 * 2 * np.pi)
    x2 = l1 * np.sin(theta1 * 2 * np.pi) + l2 * np.sin(theta2 * 2 * np.pi)
    y1 = -l1 * np.cos(theta1 * 2 * np.pi)
    y2 = -l1 * np.cos(theta1 * 2 * np.pi) -l2 * np.cos(theta2 * 2 * np.pi)
    return [x1, x2, y1, y2]

# taylor series approx

def approx(s):
    theta1 = s[0] + s[2] * deltat + s[4] / 2 * deltat**2
    theta2 = s[1] + s[3] * deltat + s[5] / 2 * deltat**2
    dt1 = s[2] + s[4] * deltat
    dt2 = s[3] + s[5] * deltat

    m1, m2, l1, l2, g = s[6], s[7], s[8], s[9], s[10]

    ddt1 = (-g *(2 * m1 + m2) *np.sin(theta1 * 2 * np.pi) - m2*g*np.sin((theta1 - 2 * theta2) * 2 * np.pi) - 2 * np.sin((theta1-theta2) * 2 * np.pi) * m2 * (dt2**2 * l2 + dt1**2 * l1 * np.cos((theta1-theta2) * 2 * np.pi))) / (l1 * (2 * m1 + m2 - m2 * np.cos((2 * theta1 - 2 * theta2) * 2 * np.pi)))

    ddt2 = (2 * np.sin((theta1 - theta2) * 2 * np.pi) * (dt1**2 * l1 * (m1 + m2) + g * (m1 + m2) * np.cos((theta1) * 2 * np.pi) + dt2**2 * l2 * m2 * np.cos((theta1-theta2) * 2 * np.pi))) / (l2 * (2 * m1 + m2 - m2 * np.cos((2 * theta1 - 2 * theta2) * 2 * np.pi)))

    return [theta1, theta2, dt1, dt2, ddt1, ddt2, m1, m2, l1, l2, g]

# updating pendulum data

def animate(history):
    for x in history:
        x1, x2, y1, y2 = get_coords(x[0], x[1], x[8], x[9])
        line.set_data([0, x1, x2], [0, y1, y2])
        plt.pause(0.02)

def generate_data(s, history = []):
    for i in t:
        history.append(s)
        s = approx(s)
    return history

# configuring figure, line parameters

fig, ax = plt.subplots(1,1, figsize=(8,8))
line, = plt.plot([], [], 'wo-', lw = 1, markersize = 5)
ax.set_xlim(-4,4)
ax.set_ylim(-4,4)
ax.set_facecolor('k')
plt.title("A Double Pendulum")

# initial parameters
m1, m2, l1, l2, g = 1, 1, 3, 1, 9.81 # (kg, kg, m, m, m/s^2)
theta1, theta2, dt1, dt2 = .25, .5, 0, 0 # (_, _, Hz, Hz)

# each theta dot dot from Lagrangian

ddt1 = (-g *(2 * m1 + m2) *np.sin(theta1 * 2 * np.pi) - m2*g*np.sin((theta1 - 2 * theta2) * 2 * np.pi) - 2 * np.sin((theta1-theta2) * 2 * np.pi) * m2 * (dt2**2 * l2 + dt1**2 * l1 * np.cos((theta1-theta2) * 2 * np.pi))) / (l1 * (2 * m1 + m2 - m2 * np.cos((2 * theta1 - 2 * theta2) * 2 * np.pi)))

ddt2 = (2 * np.sin((theta1 - theta2) * 2 * np.pi) * (dt1**2 * l1 * (m1 + m2) + g * (m1 + m2) * np.cos((theta1) * 2 * np.pi) + dt2**2 * l2 * m2 * np.cos((theta1-theta2) * 2 * np.pi))) / (l2 * (2 * m1 + m2 - m2 * np.cos((2 * theta1 - 2 * theta2) * 2 * np.pi)))

# making sliders

S = [theta1, theta2, dt1, dt2, ddt1, ddt2, m1, m2, l1, l2, g]

fig.subplots_adjust(left=0.25, bottom=0.4)

axm1 = fig.add_axes([0.25, 0.30, 0.65, 0.03])
axm2 = fig.add_axes([0.25, 0.25, 0.65, 0.03])
axl1 = fig.add_axes([0.25, 0.20, 0.65, 0.03])
axl2 = fig.add_axes([0.25, 0.15, 0.65, 0.03])
axth1 = fig.add_axes([0.25, 0.10, 0.65, 0.03])
axth2 = fig.add_axes([0.25, 0.05, 0.65, 0.03])

m1_slider = Slider(ax=axm1, label=r'$m_1$', valmin=0.001, valmax=2, valinit=m1)
m2_slider = Slider(ax=axm2, label=r'$m_2$', valmin=0.001, valmax=2, valinit=m2)
l1_slider = Slider(ax=axl1, label=r'$L_1$', valmin=0.01, valmax=2, valinit=l1)
l2_slider = Slider(ax=axl2, label=r'$L_2$', valmin=0.01, valmax=2, valinit=l2)
th1_slider = Slider(ax=axth1, label=r'$\theta_1$', valmin=0.001, valmax=1, valinit=theta1)
th2_slider = Slider(ax=axth2, label=r'$\theta_2$', valmin=0.001, valmax=1, valinit=theta2)

# make button to reset sliders and begin sim
resetax = fig.add_axes([0.09, 0.12, 0.1, 0.04])
startax = fig.add_axes([0.09, 0.07, 0.1, 0.04])

button_reset = Button(resetax, 'Reset', hovercolor='0.975')
button_start = Button(startax, 'Start', hovercolor='0.975')

def reset(event):
    m1_slider.reset()
    m2_slider.reset()
    l1_slider.reset()
    l2_slider.reset()
    th1_slider.reset()
    th2_slider.reset()
button_reset.on_clicked(reset)

def start(event):
    history = generate_data([th1_slider.val, th2_slider.val, dt1, dt2, ddt1, ddt2, m1_slider.val, m2_slider.val, l1_slider.val, l2_slider.val, g])
    animate(history)
button_start.on_clicked(start)

plt.show()


###############################################################
# ----- ----- ----- ----- Second Part ----- ----- ----- ----- #
###############################################################

# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.widgets import Slider, Button

# # each frame length

# deltat = .01
# t = np.linspace(0, 40, 1001)

# # get x,y for each mass

# def get_coords(theta1, theta2, l1, l2):
#     x1 = l1 * np.sin(theta1 * 2 * np.pi)
#     x2 = l1 * np.sin(theta1 * 2 * np.pi) + l2 * np.sin(theta2 * 2 * np.pi)
#     y1 = -l1 * np.cos(theta1 * 2 * np.pi)
#     y2 = -l1 * np.cos(theta1 * 2 * np.pi) -l2 * np.cos(theta2 * 2 * np.pi)
#     return [x1, x2, y1, y2]

# # taylor series approx

# def approx(s):
#     theta1 = s[0] + s[2] * deltat + s[4] / 2 * deltat**2
#     theta2 = s[1] + s[3] * deltat + s[5] / 2 * deltat**2
#     dt1 = s[2] + s[4] * deltat
#     dt2 = s[3] + s[5] * deltat

#     m1, m2, l1, l2, g = s[6], s[7], s[8], s[9], s[10]

#     ddt1 = (-g *(2 * m1 + m2) *np.sin(theta1 * 2 * np.pi) - m2*g*np.sin((theta1 - 2 * theta2) * 2 * np.pi) - 2 * np.sin((theta1-theta2) * 2 * np.pi) * m2 * (dt2**2 * l2 + dt1**2 * l1 * np.cos((theta1-theta2) * 2 * np.pi))) / (l1 * (2 * m1 + m2 - m2 * np.cos((2 * theta1 - 2 * theta2) * 2 * np.pi)))

#     ddt2 = (2 * np.sin((theta1 - theta2) * 2 * np.pi) * (dt1**2 * l1 * (m1 + m2) + g * (m1 + m2) * np.cos((theta1) * 2 * np.pi) + dt2**2 * l2 * m2 * np.cos((theta1-theta2) * 2 * np.pi))) / (l2 * (2 * m1 + m2 - m2 * np.cos((2 * theta1 - 2 * theta2) * 2 * np.pi)))

#     return [theta1, theta2, dt1, dt2, ddt1, ddt2, m1, m2, l1, l2, g]

# # updating pendulum data

# def animate(history):
#     theta1s = [x[0] % 2 - 1 for x in history]
#     theta2s = [x[1] % 2 - 1 for x in history]
#     for i in range(len(history)):
#         graph.set_data(theta1s[:i+1], theta2s[:i+1])

# def generate_data(s, history = []):
#     for i in t:
#         history.append(s)
#         s = approx(s)
#     return history

# # configuring figure, line parameters

# fig, ax = plt.subplots(1,1, figsize=(8,8))
# graph, = plt.plot([], [], 'ko', markersize = 1)
# ax.set_xlim(-5,5)
# ax.set_ylim(-5,5)
# plt.title("A Double Pendulum")

# # initial parameters
# m1, m2, l1, l2, g = 1, 0.5, .5, 1, 9.81 # (kg, kg, m, m, m/s^2)
# theta1, theta2, dt1, dt2 = .318, .127, 0, 0 # (_, _, Hz, Hz)

# # each theta dot dot from Lagrangian

# ddt1 = (-g *(2 * m1 + m2) *np.sin(theta1 * 2 * np.pi) - m2*g*np.sin((theta1 - 2 * theta2) * 2 * np.pi) - 2 * np.sin((theta1-theta2) * 2 * np.pi) * m2 * (dt2**2 * l2 + dt1**2 * l1 * np.cos((theta1-theta2) * 2 * np.pi))) / (l1 * (2 * m1 + m2 - m2 * np.cos((2 * theta1 - 2 * theta2) * 2 * np.pi)))

# ddt2 = (2 * np.sin((theta1 - theta2) * 2 * np.pi) * (dt1**2 * l1 * (m1 + m2) + g * (m1 + m2) * np.cos((theta1) * 2 * np.pi) + dt2**2 * l2 * m2 * np.cos((theta1-theta2) * 2 * np.pi))) / (l2 * (2 * m1 + m2 - m2 * np.cos((2 * theta1 - 2 * theta2) * 2 * np.pi)))

# # making sliders

# S = [theta1, theta2, dt1, dt2, ddt1, ddt2, m1, m2, l1, l2, g]

# fig.subplots_adjust(left=0.25, bottom=0.4)

# axm1 = fig.add_axes([0.25, 0.30, 0.65, 0.03])
# axm2 = fig.add_axes([0.25, 0.25, 0.65, 0.03])
# axl1 = fig.add_axes([0.25, 0.20, 0.65, 0.03])
# axl2 = fig.add_axes([0.25, 0.15, 0.65, 0.03])
# axth1 = fig.add_axes([0.25, 0.10, 0.65, 0.03])
# axth2 = fig.add_axes([0.25, 0.05, 0.65, 0.03])

# m1_slider = Slider(ax=axm1, label=r'$m_1$', valmin=0.001, valmax=2, valinit=m1)
# m2_slider = Slider(ax=axm2, label=r'$m_2$', valmin=0.001, valmax=2, valinit=m2)
# l1_slider = Slider(ax=axl1, label=r'$L_1$', valmin=0.01, valmax=2, valinit=l1)
# l2_slider = Slider(ax=axl2, label=r'$L_2$', valmin=0.01, valmax=2, valinit=l2)
# th1_slider = Slider(ax=axth1, label=r'$\theta_1$', valmin=0.001, valmax=1, valinit=theta1)
# th2_slider = Slider(ax=axth2, label=r'$\theta_2$', valmin=0.001, valmax=1, valinit=theta2)

# # make button to reset sliders and begin sim
# resetax = fig.add_axes([0.09, 0.12, 0.1, 0.04])
# startax = fig.add_axes([0.09, 0.07, 0.1, 0.04])

# button_reset = Button(resetax, 'Reset', hovercolor='0.975')
# button_start = Button(startax, 'Start', hovercolor='0.975')

# def reset(event):
#     m1_slider.reset()
#     m2_slider.reset()
#     l1_slider.reset()
#     l2_slider.reset()
#     th1_slider.reset()
#     th2_slider.reset()
# button_reset.on_clicked(reset)

# def start(event):
#     history = generate_data([th1_slider.val, th2_slider.val, dt1, dt2, ddt1, ddt2, m1_slider.val, m2_slider.val, l1_slider.val, l2_slider.val, g])
#     animate(history)
# button_start.on_clicked(start)

# plt.show()