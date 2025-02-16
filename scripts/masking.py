import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, RectangleSelector, Button
from astropy.io import fits
import matplotlib.colors as mcolors
from scipy.ndimage import binary_dilation

class StandaloneSegmentationEditor:
    def __init__(self, detection_hdus, seg_hdu, obj_id):
        """
        Opens a standalone segmentation editor in a separate window.
        Args:
            detection_hdus (list): List of detection image HDUs.
            seg_hdu (HDU): Segmentation map HDU.
            obj_id (int): Object ID for segmentation.
        """
        self.detection_hdus = detection_hdus
        self.seg_hdu = seg_hdu
        self.obj_id = obj_id

        # Store original segmentation for resetting
        self.original_segmentation_map = seg_hdu.data.copy()
        self.segmentation_map = self.original_segmentation_map.copy()
        self.header = seg_hdu.header  # Keep header metadata

        # ✅ Initialize `result` to prevent AttributeError
        self.result = None  

        # ✅ Create a brand new figure in a separate window
        plt.ion()  # Enable interactive mode
        self.fig = plt.figure(figsize=(15, 5))  # New figure, forces separate window
        self.axes = [self.fig.add_subplot(1, 4, i+1) for i in range(4)]  # 4 Subplots

        self.canvas = self.fig.canvas

        # ✅ Automatically maximize the window
        try:
            self.fig_manager = plt.get_current_fig_manager()
            self.fig_manager.window.showMaximized()
        except AttributeError:
            print("Warning: Unable to maximize window on some systems.")

        self.load_images()  # Load images
        self.init_widgets()  # Setup interactive tools

        plt.show(block=True)  # ✅ Open in a truly separate window

    def load_images(self):
        """Load and display images with segmentation overlays using the same colormap but grayscale outside seg."""
        self.image_plots = []
        seg_mask = self.segmentation_map == self.obj_id

        for i in range(3):
            img = self.detection_hdus[i].data

            # ✅ Normalize data to maintain same color scale
            norm = mcolors.Normalize(vmin=np.nanmin(img), vmax=np.nanmax(img))

            # ✅ Apply plasma_r colormap to both, but make outside grayscale
            seg_only = np.where(seg_mask, img, np.nan)
            non_seg_only = np.where(~seg_mask, img, np.nan)

            # ✅ Convert grayscale using `to_rgba()` with plasma_r but desaturate outside
            cmap = plt.cm.plasma_r
            rgba_img = cmap(norm(img))  # Convert to RGBA
            rgba_img[~seg_mask, :3] = np.mean(rgba_img[~seg_mask, :3], axis=1, keepdims=True)  # Convert outside to grayscale

            non_seg_plot = self.axes[i].imshow(rgba_img, origin="lower")
            seg_plot = self.axes[i].imshow(seg_only, cmap="plasma_r", origin="lower", norm=norm)

            self.axes[i].set_title(f"Detection Image {i+1}")
            self.image_plots.append((non_seg_plot, seg_plot))

        # ✅ Show segmentation map
        self.seg_im = self.axes[3].imshow(self.segmentation_map, cmap="jet", origin="lower")

        self.axes[3].set_title("Segmentation Map with Pixel Boundaries")
        self.canvas.draw()

    def update_images(self):
        """Update segmentation overlay with the same plasma_r colormap but grayscale outside segmentation."""
        seg_mask = self.segmentation_map == self.obj_id


        # ✅ Clear previous contours before adding new ones
        for ax in self.axes:
            for coll in ax.collections:
                coll.remove()

        for i in range(3):
            img = self.detection_hdus[i].data
            norm = mcolors.Normalize(vmin=np.nanmin(img), vmax=np.nanmax(img))

            # ✅ Convert plasma_r to grayscale for non-segmented areas
            rgba_img = plt.cm.plasma_r(norm(img))
            rgba_img[~seg_mask, :3] = np.mean(rgba_img[~seg_mask, :3], axis=1, keepdims=True)  # Convert outside to grayscale

            self.image_plots[i][0].set_data(rgba_img)  # ✅ Update grayscale non-seg region
            self.image_plots[i][1].set_data(np.where(seg_mask, img, np.nan))  # ✅ Update segmented region
        self.seg_im.set_data(self.segmentation_map)  # ✅ Update segmentation map

        self.canvas.draw_idle()  # ✅ Efficient update


    def modify_segmentation(self, x1, x2, y1, y2, remove_only=False, add_only=False):
        """Modify pixels in the selected region based on mouse action."""
        x1, x2 = int(min(x1, x2)), int(max(x1, x2))
        y1, y2 = int(min(y1, y2)), int(max(y1, y2))

        if remove_only:
            self.segmentation_map[y1:y2, x1:x2][self.segmentation_map[y1:y2, x1:x2] == self.obj_id] = 0
        elif add_only:
            self.segmentation_map[y1:y2, x1:x2][self.segmentation_map[y1:y2, x1:x2] == 0] = self.obj_id
        else:
            mask = self.segmentation_map[y1:y2, x1:x2] == self.obj_id
            self.segmentation_map[y1:y2, x1:x2] = np.where(mask, 0, self.obj_id)

        self.update_images()

    def on_click(self, event):
        """Handles mouse clicks for modifying segmentation."""
        if event.inaxes is None:
            return

        x, y = int(event.xdata), int(event.ydata)
        if event.button == 1:  # Left click: Remove
            self.modify_segmentation(x, x + 1, y, y + 1, remove_only=True)
        elif event.button == 3:  # Right click: Add
            self.modify_segmentation(x, x + 1, y, y + 1, add_only=True)

    def on_select(self, eclick, erelease):
        """Handles rectangle selection."""
        if eclick.inaxes is None:
            return

        x1, y1 = int(eclick.xdata), int(eclick.ydata)
        x2, y2 = int(erelease.xdata), int(erelease.ydata)

        if eclick.button == 1:  # Left button drag: Remove segmentation
            self.modify_segmentation(x1, x2, y1, y2, remove_only=True)
        elif eclick.button == 3:  # Right button drag: Add segmentation
            self.modify_segmentation(x1, x2, y1, y2, add_only=True)

    def reset_segmentation(self, event):
        """Reset segmentation map to original."""
        self.segmentation_map = self.original_segmentation_map.copy()
        self.update_images()

    def finish(self, event):
        """Close figure and store segmentation result."""
        self.result = fits.ImageHDU(data=self.segmentation_map, header=self.header)  # ✅ Store result before closing
        plt.close(self.fig)

    def close_event(self, event):
        """Ensure result is stored even if window is closed manually."""
        if self.result is None:  # ✅ If "Finish" was not clicked, store the result
            self.result = fits.ImageHDU(data=self.segmentation_map, header=self.header)

    def init_widgets(self):
        """Initialize interactive tools (cursor, rectangle selector, buttons)."""
        self.cursors = []
        self.rect_selectors = []

        for ax in self.axes:
            # ✅ Apply crosshair cursor to all panels
            cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
            self.cursors.append(cursor)

            # ✅ Apply rectangle selection tool to all panels
            selector = RectangleSelector(
                ax, self.on_select,
                useblit=True, interactive=True,
                props=dict(facecolor='red', edgecolor='black', alpha=0.3, fill=True)
            )
            selector.set_active(True)
            self.rect_selectors.append(selector)

        # Buttons for reset & finish
        reset_ax = self.fig.add_axes([0.75, 0.02, 0.1, 0.05])
        finish_ax = self.fig.add_axes([0.85, 0.02, 0.1, 0.05])

        self.reset_button = Button(reset_ax, 'Reset')
        self.finish_button = Button(finish_ax, 'Finish')

        self.reset_button.on_clicked(self.reset_segmentation)
        self.finish_button.on_clicked(self.finish)

        self.canvas.mpl_connect("button_press_event", self.on_click)  # Mouse click event
        self.fig.canvas.mpl_connect("close_event", self.close_event)  # ✅ Ensure result is stored if closed manually

    def get_result(self):
        """Return final segmentation as ImageHDU."""
        return self.result


def run_segmentation_editor(detection_hdus, seg_hdu, obj_id):
    """Run the segmentation editor in a separate window."""
    editor = StandaloneSegmentationEditor(detection_hdus, seg_hdu, obj_id)
    return editor.get_result()
