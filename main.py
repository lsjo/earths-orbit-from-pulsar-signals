from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

def parabola(x1, y1, x2, y2, x3, y3, plot=False):
    denom = (x1-x2) * (x1-x3) * (x2-x3)
    A = (x3 * (y2-y1) + x2 * (y1-y3) + x1 * (y3-y2)) / denom
    B = (x3*x3 * (y1-y2) + x2*x2 * (y3-y1) + x1*x1 * (y2-y3)) / denom
    C = (x2 * x3 * (x2-x3) * y1+x3 * x1 * (x3-x1) * y2+x1 * x2 * (x1-x2) * y3)/ denom
    if abs(x3-x2) < 1 or abs(x2-x1) < 1:
        print(x1, x2, x3)
    if plot == True:
        xList = []
        yList = []
        for n in range(int(x1), int(x3), 1):
            xList.append(n)
            yList.append(A*n*n + B*n + C)
            plt.plot(xList, yList, c="k")
        return A, B, C
    else:
        return A,B,C

def distance(x1, y1, x2, y2):
    return np.sqrt((x1-x2)** 2 + (y1-y2)**2)

def delta(df, column, offset):
    t0 = df.iloc[0][column]
    index = 0
    valsSum = 0
    averages = []
    while (t0-df.iloc[index][column] < offset) and (index < len(df[column]) and (1==0)):
        valsSum = valsSum + (df.iloc[index][column])
    averages.append(valsSum)
    return df, averages
        

def interpolate(df, c1, c1Value, c2, method="linear"):
    df.sort_values(by=[c1])
    df.reset_index(drop=True)
    
    for n in df[c1]:
        value = n
        if n > c1Value:
            break
    valueIndex = df.isin([value]).any(axis=1).idxmax()
    
    x1 = df.iloc[valueIndex-1][c1]
    x2 = df.iloc[valueIndex][c1]
    x3 = df.iloc[valueIndex+1][c1]
    y1 = df.iloc[valueIndex-1][c2]
    y2 = df.iloc[valueIndex][c2]
    y3 = df.iloc[valueIndex+1][c2]

    if method == "linear":
        k = (y2-y1)/(x2-x1)
    
        interval = (c1Value - x1)/(x2 - x1)
        c2Value = y1 + interval*(y2-y1)
    else:
        a, b, c = parabola(x1, y1, x2, y2, x3, y3, plot=True)
        c2Value = (a * (c1Value ** 2)) + (b * c1Value) + c

    return c2Value

print("-----------------")
print("PULSAR J1614-2230")
print("-----------------")
print()
J1614 = pd.read_csv("/data/J1614-2230.residuals.csv")
J1614s = pd.read_csv("/data/J1614-2230.residuals.SIMPLIFIED.csv")
print(J1614.head())
print()
print(J1614["DELAY"].describe())
print()

print("-----------------")
print("PULSAR J2222-0137")
print("-----------------")
print()
J2222 = pd.read_csv("/data/J2222-0137.residuals.csv")
J2222s = pd.read_csv("/data/J2222-0137.residuals.SIMPLIFIED.csv")
print(J2222.head())
print()
print(J2222["DELAY"].describe())
print()


J1614["DELAY"] = J1614["DELAY"] * np.cos(-1.257 * (np.pi/180))
J2222["DELAY"] = J2222["DELAY"] * np.cos(7.977 * (np.pi/180))

J1614s["DELAY"] = J1614s["DELAY"] * np.cos(-1.257 * (np.pi/180))
J2222s["DELAY"] = J2222s["DELAY"] * np.cos(7.977 * (np.pi/180))


J1614["DELAY A"] = J1614["DELAY"]-((J1614["DELAY"].max()+J1614["DELAY"].min())/2)
J2222["DELAY A"] = J2222["DELAY"]-((J2222["DELAY"].max()+J2222["DELAY"].min())/2)

J1614s["DELAY A"] = J1614s["DELAY"]-((J1614s["DELAY"].max()+J1614s["DELAY"].min())/2)
J2222s["DELAY A"] = J2222s["DELAY"]-((J2222s["DELAY"].max()+J2222s["DELAY"].min())/2)

plt.scatter(x=J1614["DAY"], y=J1614["DELAY"])
plt.title("Delay for PSR J1614-2230 vs Time")
plt.xlabel("Time of Pulse (day)")
plt.ylabel("Delay of Pulse (ms)")
plt.axhline(J1614["DELAY"].median(), label="median", c="b")
plt.axhline(J1614["DELAY"].mean(), label="mean", c="r")
plt.axhline((J1614["DELAY"].max()+J1614["DELAY"].min())/2, label="(max-min)/2", c="g")
plt.legend()
plt.show()

plt.scatter(x=J2222["DAY"], y=J2222["DELAY"])
plt.title("Delay for PSR J2222-0137 vs Time")
plt.xlabel("Time of Pulse (day)")
plt.ylabel("Delay of Pulse (ms)")
plt.axhline(J2222["DELAY"].median(), label="median", c="b")
plt.axhline(J2222["DELAY"].mean(), label="mean", c="r")
plt.axhline((J2222["DELAY"].max()+J2222["DELAY"].min())/2, label="(max-min)/2", c="g")
plt.legend()
plt.show()

# plt.plot(J1614["DAY"], J1614["DELAY"]-((J1614["DELAY"].max()+J1614["DELAY"].min())/2), c="b", linestyle=":")
# plt.plot(J2222["DAY"], J2222["DELAY"]-((J2222["DELAY"].max()+J2222["DELAY"].min())/2), c="r", linestyle=":")
plt.scatter(J1614["DAY"], J1614["DELAY A"], c="b", label="J1614")
plt.scatter(J2222["DAY"], J2222["DELAY A"], c="r", label="J2222")
plt.title("Delay vs Time")
plt.xlabel("Time of Pulse (day)")
plt.ylabel("Delay of Pulse (ms)")
plt.legend(loc="lower left")
plt.show()

J1614 = J1614[J1614["DAY"] >= 58800]
J1614 = J1614[J1614["DAY"] <= 60200]
J1614 = J1614.reset_index()

J2222 = J2222[J2222["DAY"] >= 58800]
J2222 = J2222[J2222["DAY"] <= 60200]
J2222 = J2222.reset_index()          

interpVals = []

shortenedTimes = J1614s.iloc[0:-2]["DAY"]

for n in shortenedTimes:
    interpVals.append(interpolate(J2222s, "DAY", n, "DELAY", method="linear"))

for n in range(2):
    interpVals.append(np.nan)

J1614s["DELAY (INTERP)"] = interpVals

J1614s = J1614s.drop(J1614s[J1614s['DELAY (INTERP)'] > 1000].index)
J1614s = J1614s.drop(J1614s[J1614s['DELAY (INTERP)'] < -1500].index)

J1614s, a = delta(J1614s,"DELAY", 1)

plt.scatter(J1614s["DAY"], J1614s["DELAY"], c="b", label="J1614")
plt.scatter(J1614s.iloc[0:-3]["DAY"], J1614s.iloc[0:-3]["DELAY (INTERP)"], c="g", label="J2222 Interpolated")
plt.scatter(J2222s["DAY"], J2222s["DELAY"], c="r", label="J2222 Raw")
plt.legend()
plt.xlabel("Day")
plt.ylabel("Delay")
plt.show()

plt.scatter(J1614s[0:-2]["DELAY"], J1614s[0:-2]["DELAY (INTERP)"], c="b", label="Earth's Orbit")
plt.legend()
plt.xlabel("J1614 Delay (ms)")
plt.xscale("linear")
plt.ylabel("J2222 Delay (ms)")
plt.axis("equal")
plt.title("Earth's orbit from J1614 and J2222 Timing Delays")
plt.show()

yearCategories = [0, 58954.95232, 59319.95232, 59684.95232, 60049.95232]
colours = [None, "#00648f", "#a234eb", "#eb3434", "#0fa84c"]

for n in yearCategories:
    if n == 0:
        continue
    index = yearCategories.index(n)
    Label = "Year " + str(index)
    plt.scatter(J1614s[(J1614s["DAY"] <= n) & (J1614s["DAY"] >= yearCategories[index-1])][0:-2]["DELAY"],
                J1614s[(J1614s["DAY"] <= n) & (J1614s["DAY"] >= yearCategories[index-1])][0:-2]["DELAY (INTERP)"],
                c=colours[index], label=Label)

plt.legend(loc = "lower left")
plt.xlabel("J1614 Delay (ms)")
plt.ylabel("J2222 Delay (ms)")
plt.axis("equal")
plt.title("Earth's orbit from J1614 and J2222 Timing Delays")
plt.show()

J1614s = J1614s.sort_values(by=["DELAY"])
J1614sMin = J1614s.iloc[0]["DELAY"]
J1614sMax = J1614s.iloc[-1]["DELAY"]
diameterX = J1614sMax - J1614sMin
radiusX = diameterX / 2

J2222s = J2222s.sort_values(by=["DELAY"])
J2222sMin = J2222s.iloc[0]["DELAY"]
J2222sMax = J2222s.iloc[-1]["DELAY"]
diameterY = J2222sMax - J2222sMin
radiusY = diameterY / 2

J1614s = J1614s.sort_values(by=["DAY"])

speedList = []
distanceList = []
timeList = []

for n in range(len(J1614s["DAY"])-2):
    x1 = J1614s.iloc[n]["DELAY"]
    y1 = J1614s.iloc[n]["DELAY (INTERP)"]
    x2 = J1614s.iloc[n+1]["DELAY"]
    y2 = J1614s.iloc[n+1]["DELAY (INTERP)"]
    Distance = distance(x1, y1, x2, y2)
    time = J1614s.iloc[n+1]["DAY"] - J1614s.iloc[n]["DAY"]
    distanceList.append(Distance)
    timeList.append(time)
    speedList.append(Distance/time)

for n in range(2):
    speedList.append(np.nan)
    timeList.append(np.nan)
    distanceList.append(np.nan)

J1614s["SPEED"]= speedList
J1614s["SPEED"] = J1614s["SPEED"] * 2.9979 * 10**8 / 86400 / 1000

J1614sd = J1614s[J1614s["DAY"] > 58700]

plt.plot(J1614sd["DAY"], J1614sd["SPEED"])
plt.xlabel("Days")
plt.ylabel("Velocity (km/s)")
plt.title("Velocity vs Time for Earth's Orbit")
plt.show()

X = J1614sd["SPEED"]** 4 * J1614sd["DELAY"] / 500 * 10
Y = J1614sd["SPEED"]** 4 * J1614sd["DELAY (INTERP)"] / 500 * 10

plt.scatter(X, Y)
plt.title("Speed of Earth in Orbit")
plt.xlabel("Earth's Speed^4 * J1614 Delay")
plt.ylabel("Earth's Speed^4 * J2222 Delay")
plt.axis("equal")
plt.show()
