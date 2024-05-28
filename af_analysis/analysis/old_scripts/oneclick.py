import matplotlib.pyplot as plt

def create_annot(ax, sc, labels, ind):
    pos = sc.get_offsets()[ind["ind"][0]]
    annot = ax.annotate("", xy=pos, xytext=(20,20), textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    text = ', '.join([labels[n] for n in ind['ind']])
    annot.set_text(text)
    return annot

def onclick(event, ax, sc, annot_list, labels):
    if event.inaxes == ax:
        cont, ind = sc.contains(event)
        if cont:
            clicked_pos = sc.get_offsets()[ind["ind"][0]]

            # Check if an annotation for this point already exists
            for annot in annot_list:
                # Compare each element of the positions
                if all(annot.xy == clicked_pos):
                    annot.set_visible(False)  # Hide the annotation
                    annot_list.remove(annot)  # Remove it from the list
                    plt.draw()
                    return

            # If not, create a new annotation
            annot = create_annot(ax, sc, labels, ind)
            annot.set_visible(True)
            annot_list.append(annot)
            plt.draw()

def plot_scatter(X, Y, x_label, y_label, custom_labels):
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 18,
           }

    fig, ax = plt.subplots(figsize=(5, 5))  # Set the figure size
    sc = ax.scatter(X, Y, s=2)  # Create a scatter plot
    ax.set_xlabel(x_label, fontdict=font)  # Set the x-axis label
    ax.set_ylabel(y_label, fontdict=font)  # Set the y-axis label
    ax.tick_params(axis='both', labelsize=font['size'])  # Set tick label sizes
    plt.tight_layout()  # Adjust the layout to prevent clipping

    annot_list = []
    fig.canvas.mpl_connect("button_press_event", lambda event: onclick(event, ax, sc, annot_list, custom_labels))

    plt.show()

# Sample data
X = [1, 2, 3, 4, 5]
Y = [2, 3, 5, 7, 11]
custom_labels = ['label1', 'label2', 'label3', 'label4', 'label5']  # Custom labels for each data point
plot_scatter(X, Y, "X-axis", "Y-axis", custom_labels)
