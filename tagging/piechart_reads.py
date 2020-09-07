

# TODO: calculate read count of each part
# tide_ value read count from sam files
p1 = 20634-(11125 + 4510 + 1337)
p2 = 11125 - (4510 + 1337)
p3 = 4510 - 1337
p4 = 1337

vals = np.array([[275., 1223., 2089., 9666], [p4, p3, p2, p1], [16779, 295, 0 , 0]])
vals.sum(axis=1)

x = np.linspace(0.0, 1.0, 12)

fig, ax = plt.subplots(subplot_kw=dict(aspect="equal"))
size = 0.3
vals = np.array([[275., 1223., 2089., 9666], [p4, p3, p2, p1], [16779, 295, 0 , 0]])

cmap = plt.get_cmap("tab20c") # only integer.
outer_colors = cmap(np.array([0,8,16])) # blue: with backbone # green without backbone
inner_colors = cmap(np.array([0, 1, 2,3,  8,9, 10, 11, 17,18, 18, 18]))
# outer
wedges_inner, texts= ax.pie(vals.sum(axis=1), radius=1-size, colors=outer_colors,
       wedgeprops=dict(width=size, edgecolor='w'))


# inner circle

float_text = ["20+ repeats", "10-20 repeats", "6-10 repeats", "1-5 repeats",
              "20+ repeats", "10-20 repeats", "6-10 repeats", "1-5 repeats",
              "very short inserts (<50bp)","only backbones","",""]

wedges, texts = ax.pie(vals.flatten(), radius=1, colors=inner_colors,
       wedgeprops=dict(width=size, edgecolor='w'))

data = vals.flatten()

bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"),
          bbox=bbox_props, zorder=0, va="center")

for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax.annotate(float_text[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
                horizontalalignment=horizontalalignment, **kw)

ax.set(aspect="equal", title='Cyclomics read information distribution')
ax.set_title("P260(DER4458)")
ax.legend(wedges_inner, ["read with backbone", "read without backbone", "only backbone (not informative)"],
          title="Total 69235 reads",
          loc="center left",
          bbox_to_anchor=(1, -1, 0.5, 1))

plt.show()

