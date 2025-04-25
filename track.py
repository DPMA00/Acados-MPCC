
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import splprep, splev
import pandas as pd
from scipy.interpolate import CubicSpline



class Track:
    def __init__(self, data='track_control_pts.csv'):
        self.df = pd.read_csv(data, delimiter=',')
        xy = self.df[['x', 'y']].values
        x = xy[:, 0]
        y = xy[:, 1]
        tck, u = splprep(xy.T, u=None, s=0.0, per=1)


        u_new = np.linspace(u.min(), u.max(), 1000)
        self.x_new, self.y_new = splev(u_new, tck, der=0)

        dx = np.diff(self.x_new)
        dy = np.diff(self.y_new)

        s = np.cumsum(np.sqrt(dx**2+dy**2))
        s = np.insert(s,0,0)
        self.L = s[-1]


        self.cs_x = CubicSpline(s, self.x_new, bc_type='periodic')
        self.cs_y = CubicSpline(s, self.y_new, bc_type='periodic')


    def create_grid(self, density=5000):

        theta_p = np.linspace(0,self.L, density)

        x_p = self.cs_x(theta_p)
        y_p = self.cs_y(theta_p)
        
        return x_p, y_p, theta_p




    def plot_track(self):
        plt.figure(figsize=(8, 6))
        plt.plot(self.x_new, self.y_new, label='Racetrack Path')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.grid(True)
        plt.axis('equal')
        plt.legend()
        plt.show()

    def get_track_plot_params(self):
        return self.x_new, self.y_new



