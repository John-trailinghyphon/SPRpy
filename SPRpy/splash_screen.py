# Start menu and splash screen of SPRpy: an SPR curve (calculated with SPRpy's own Fresnel model), buttons to
# open a measurement, reopen a session or convert .spr2 files, and the loading status while SPRpy starts.
# The start menu is optional (show_start_menu in config.toml) and never stops SPRpy from starting: any problem with the
# display falls back to the console prompts.

import glob
import os
import time
import traceback
import numpy as np

# Console banner, assembled letter by letter ("big" figlet font)
_LETTERS = {
    'S': ['  _____ ', ' / ____|', '| (___  ', r' \___ \ ', ' ____) |', r'|_____/ ', '        ', '        '],
    'P': [' _____  ', r'|  __ \ ', '| |__) |', '|  ___/ ', '| |     ', '|_|     ', '        ', '        '],
    'R': [' _____  ', r'|  __ \ ', '| |__) |', '|  _  / ', r'| | \ \ ', '|_|  \\_\\', '        ', '        '],
    'p': ['       ', '       ', ' _ __  ', r"| '_ \ ", '| |_) |', '| .__/ ', '| |    ', '|_|    '],
    'y': ['       ', '       ', ' _   _ ', '| | | |', '| |_| |', r' \__, |', '  __/ |', ' |___/ '],
}

WIDTH = 640
SPLASH_HEIGHT, MENU_HEIGHT = 360, 470
TOP_COLOR, BOTTOM_COLOR = (31, 52, 79), (68, 110, 155)
CURVE_COLOR, ACCENT_COLOR, TEXT_COLOR, MUTED_COLOR = '#9fd3ff', '#ffb347', '#ffffff', '#c9d6e3'
BUTTON_COLOR, BUTTON_ACTIVE_COLOR, BUTTON_TEXT_COLOR = '#ffffff', '#dbe8f5', '#2f4f73'


def banner(version, word='SPRpy'):
    # Join the letters, letting neighbouring letters share their edge column where it does not clash (like figlet)
    rows = list(_LETTERS[word[0]])
    for letter in word[1:]:
        glyph = _LETTERS[letter]
        width = max(len(row) for row in rows)
        rows = [row.ljust(width) for row in rows]
        if all(row[-1] == ' ' or g[0] == ' ' or row[-1] == g[0] for row, g in zip(rows, glyph)):
            rows = [row[:-1] + (g[0] if row[-1] == ' ' else row[-1]) + g[1:] for row, g in zip(rows, glyph)]
        else:
            rows = [row + g for row, g in zip(rows, glyph)]
    rows = [' ' + row.rstrip() for row in rows]
    rows[-1] = rows[-1].ljust(max(len(row) for row in rows) + 3) + 'v' + version
    return '\n' + '\n'.join(rows) + '\n Multi-parameter surface plasmon resonance analysis\n'


def print_banner(version):
    print(banner(version))


def _gold_spectrum(sensor_values, angles):
    """
    Reflectivity of a Cr/Au sensor in water at 670 nm with the default sensor values of config.toml.
    """
    from fresnel_transfer_matrix import fresnel_calculation
    d = np.array([np.nan, float(sensor_values['d_nm'][1]), float(sensor_values['d_nm'][2]), np.nan])
    n = np.array([float(sensor_values['n_670'][0]), float(sensor_values['n_670'][1]), float(sensor_values['n_670'][2]), 1.333])
    k = np.array([float(sensor_values['k_670'][0]), float(sensor_values['k_670'][1]), float(sensor_values['k_670'][2]), 0.0])
    return fresnel_calculation(angles=angles, wavelength=670, layer_thicknesses=d, n_re=n, n_im=k)


SUPERSAMPLING = 3  # The graphics are drawn at this many times the resolution and scaled down, for smooth (antialiased) edges


def _enable_dpi_awareness():
    """
    Tells Windows that SPRpy handles display scaling itself, so that windows are drawn at the full resolution of the
    screen instead of being drawn small and stretched (which makes text and graphics blurry at e.g. 125 % or 150 %).
    """
    try:
        import ctypes
        try:
            ctypes.windll.shcore.SetProcessDpiAwareness(1)  # System DPI aware (Windows 8.1 and later)
        except Exception:
            ctypes.windll.user32.SetProcessDPIAware()  # Older Windows versions
    except Exception:
        pass  # Not Windows


def _dashed_line(draw, x, y0, y1, dash, gap, fill, width):
    """
    Vertical dashed line from y0 to y1 (Pillow has no dashed lines).
    """
    y = y0
    while y < y1:
        draw.line([(x, y), (x, min(y + dash, y1))], fill=fill, width=width)
        y += dash + gap


def _logo_curve(svg_path):
    """
    The SPR curve of the SPRpy logo from icon.svg: points of its stroked path (viewBox units), the size of the viewBox,
    the stroke color and the stroke width. The path only uses relative move/line commands (m, l, h).
    """
    import re
    text = open(svg_path, encoding='utf-8').read()
    view_box = [float(v) for v in re.search(r'viewBox="([^"]+)"', text).group(1).split()]
    translate = re.search(r'transform="translate\(([-\d.e]+),([-\d.e]+)\)"', text)
    dx, dy = (float(translate.group(1)), float(translate.group(2))) if translate else (0.0, 0.0)
    path = max(re.finditer(r'<path\s+d="([^"]+)"(.*?)/>', text, re.S), key=lambda match: len(match.group(1)))
    color = re.search(r'stroke:(#[0-9a-fA-F]{6})', path.group(2)).group(1)
    stroke_width = float(re.search(r'stroke-width:([\d.]+)', path.group(2)).group(1))
    tokens = re.findall(r'[MLHmlh]|[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', path.group(1))
    points, x, y, command, index = [], 0.0, 0.0, 'm', 0
    while index < len(tokens):
        if tokens[index] in 'MLHmlh':
            command = tokens[index]
            index += 1
            continue
        if command in 'hH':
            value = float(tokens[index])
            x = x + value if command == 'h' else value
            index += 1
        else:
            px, py = float(tokens[index]), float(tokens[index + 1])
            x, y = (x + px, y + py) if command in 'ml' else (px, py)
            index += 2
            if command == 'm':
                command = 'l'  # Pairs after a move are lines
            elif command == 'M':
                command = 'L'
        points.append((x + dx, y + dy))
    return points, view_box[2], _hex_to_rgb(color), stroke_width


def _hex_to_rgb(color):
    return tuple(int(color[i:i + 2], 16) for i in (1, 3, 5))


class SplashScreen:
    """
    Borderless window shown while SPRpy starts. Usage:

        splash = SplashScreen(version, sensor_values)
        splash.animate()           # Shows the window for a moment before loading continues
        splash.status('Loading')   # Updates the status line
        splash.close()             # Fades out and closes

    The window is scaled to the display scaling of Windows, and the graphics (background, logo, SPR curve) are drawn
    as one antialiased image with Pillow, so that they are sharp at any scaling. Without Pillow they are drawn with the
    canvas instead.
    """

    def __init__(self, version, sensor_values=None, height=SPLASH_HEIGHT, topmost=True):
        _enable_dpi_awareness()
        import tkinter
        self.tk = tkinter
        self.height = height
        self.root = tkinter.Tk()
        self.root.title('SPRpy')
        self.root.overrideredirect(True)
        self.root.attributes('-topmost', topmost)
        self.scale = max(1.0, self.root.winfo_fpixels('1i') / 96.0)  # Display scaling of Windows (1.25 at 125 %)
        self.width_px, self.height_px = self.px(WIDTH), self.px(height)
        x = (self.root.winfo_screenwidth() - self.width_px) // 2
        y = (self.root.winfo_screenheight() - self.height_px) // 2
        self.root.geometry('{w}x{h}+{x}+{y}'.format(w=self.width_px, h=self.height_px, x=x, y=y))
        self.canvas = tkinter.Canvas(self.root, width=self.width_px, height=self.height_px, highlightthickness=0, bd=0)
        self.canvas.pack()
        self.closed = False
        self._animation_timer = None

        # SPR curve of a gold sensor (Fresnel model)
        self.plot = (60, 140, WIDTH - 40, 290)  # left, top, right, bottom (in unscaled units)
        self.angles = np.linspace(58, 78, 600)
        try:
            self.reflectivity = _gold_spectrum(sensor_values, self.angles) if sensor_values else None
        except Exception:
            self.reflectivity = None
        if self.reflectivity is None or not np.all(np.isfinite(self.reflectivity)):
            # Generic SPR-like curve if the sensor values can not be used
            self.reflectivity = 0.75 - 0.65 * np.exp(-((self.angles - 69) / 1.6) ** 2) + 0.15 * (self.angles > 61)
        self.SPR_index = int(np.argmin(self.reflectivity))
        self.TIR_index = int(np.argmax(np.diff(self.reflectivity[:self.SPR_index]))) if self.SPR_index > 2 else 0

        # Graphics: one antialiased image, or canvas items if Pillow is not available
        self._icon_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static', 'images', 'icon.png')
        self._logo_box = (32, 30, 32 + 88, 30 + 88)  # White circle behind the logo
        try:
            self._image = self._render_graphics()
            self.canvas.create_image(0, 0, image=self._image, anchor='nw', tags='background')
            title_x = self._logo_box[2] + 18
        except Exception:
            title_x = self._draw_graphics_on_canvas()

        # Text (drawn by the canvas, sharp at any scaling)
        self.text(title_x, 40, 'SPRpy', size=46, weight='bold', anchor='nw')
        self.text(title_x + 4, 106, 'Multi-parameter surface plasmon resonance analysis', size=15, color=MUTED_COLOR, anchor='nw')
        self.text(WIDTH - 46, 30, 'v' + version, size=15, weight='bold', color=MUTED_COLOR, anchor='ne')
        left, top, right, bottom = self.plot
        self.text((left + right) / 2, bottom + 14, 'Incident angle', size=12, color=MUTED_COLOR)
        self.text(left - 12, (top + bottom) / 2, 'R', size=14, slant='italic', color=MUTED_COLOR)
        x_SPR, _ = self._point(self.SPR_index)
        self.text(x_SPR + 8, top + 4, 'θ SPR', size=14, weight='bold', color=ACCENT_COLOR, anchor='nw')
        if self.TIR_index > 0:
            x_TIR, _ = self._point(self.TIR_index)
            self.text(x_TIR - 6, top + 4, 'θ TIR', size=14, color=MUTED_COLOR, anchor='ne')

        self.status_item = self.text(28, height - 26, 'Starting...', size=14, color=MUTED_COLOR, anchor='w')
        self.root.update()

    # Scaling and drawing helpers -------------------------------------------------------------------------------------
    def px(self, value):
        """
        Size or position in screen pixels from the unscaled layout units (pixels at 100 % display scaling).
        """
        return int(round(value * self.scale))

    def font(self, size, weight='normal', slant='roman'):
        return ('Segoe UI', -self.px(size), weight, slant)

    def text(self, x, y, text, size=14, weight='normal', slant='roman', color=TEXT_COLOR, anchor='center', tags=None):
        return self.canvas.create_text(self.px(x), self.px(y), text=text, anchor=anchor, fill=color,
                                       font=self.font(size, weight, slant), tags=tags)

    def _point(self, index):
        """
        Position (unscaled units) of a point of the SPR curve in the plot area.
        """
        left, top, right, bottom = self.plot
        r_max = max(1.0, float(np.max(self.reflectivity)))
        x = left + (self.angles[index] - self.angles[0]) / (self.angles[-1] - self.angles[0]) * (right - left)
        y = bottom - self.reflectivity[index] / r_max * (bottom - top)
        return x, y

    def _render_graphics(self):
        """
        Background gradient, logo, axes and SPR curve as one image, drawn with Pillow at SUPERSAMPLING times the
        resolution and scaled down for smooth edges.
        """
        from PIL import Image, ImageDraw, ImageTk
        f = self.scale * SUPERSAMPLING
        width, height = int(round(WIDTH * f)), int(round(self.height * f))

        # Vertical gradient
        fraction = np.linspace(0, 1, height)[:, np.newaxis]
        rows = np.array(TOP_COLOR) + (np.array(BOTTOM_COLOR) - np.array(TOP_COLOR)) * fraction
        image = Image.fromarray(np.repeat(rows[:, np.newaxis, :], width, axis=1).astype(np.uint8), 'RGB')
        draw = ImageDraw.Draw(image)

        # White circle with the logo curve (from icon.svg, drawn at full resolution)
        draw.ellipse([v * f for v in self._logo_box], fill=(255, 255, 255))
        logo_drawn = False
        try:
            logo_points, logo_size, logo_color, logo_stroke = _logo_curve(self._icon_path[:-4] + '.svg')
            size = 0.7 * (self._logo_box[2] - self._logo_box[0])  # Logo square inside the circle (unscaled units)
            x0 = (self._logo_box[0] + self._logo_box[2] - size) / 2
            y0 = (self._logo_box[1] + self._logo_box[3] - size) / 2
            draw.line([((x0 + px / logo_size * size) * f, (y0 + py / logo_size * size) * f) for px, py in logo_points],
                      fill=logo_color, width=max(1, int(round(logo_stroke / logo_size * size * f))), joint='curve')
            logo_drawn = True
        except Exception:
            pass

        # Axes, TIR and SPR angle markers and the SPR curve
        left, top, right, bottom = self.plot
        muted, accent, curve = _hex_to_rgb(MUTED_COLOR), _hex_to_rgb(ACCENT_COLOR), _hex_to_rgb(CURVE_COLOR)
        line = max(1, int(round(f)))
        draw.line([(left * f, bottom * f), (right * f, bottom * f)], fill=muted, width=line)
        draw.line([(left * f, top * f), (left * f, bottom * f)], fill=muted, width=line)
        if self.TIR_index > 0:
            x_TIR, _ = self._point(self.TIR_index)
            _dashed_line(draw, x_TIR * f, top * f, bottom * f, 2 * f, 4 * f, muted, line)
        x_SPR, y_SPR = self._point(self.SPR_index)
        _dashed_line(draw, x_SPR * f, top * f, bottom * f, 4 * f, 3 * f, accent, 2 * line)
        points = [(x * f, y * f) for x, y in (self._point(i) for i in range(len(self.angles)))]
        draw.line(points, fill=curve, width=int(round(3 * f)), joint='curve')
        radius = 5 * f
        draw.ellipse([x_SPR * f - radius, y_SPR * f - radius, x_SPR * f + radius, y_SPR * f + radius], fill=accent)

        image = image.resize((self.width_px, self.height_px), Image.LANCZOS)

        # Logo from icon.png if it could not be drawn from icon.svg
        try:
            if logo_drawn:
                raise StopIteration
            logo = Image.open(self._icon_path).convert('RGBA')
            size = self.px(72)
            logo.thumbnail((size, size), Image.LANCZOS)
            center_x = self.px((self._logo_box[0] + self._logo_box[2]) / 2)
            center_y = self.px((self._logo_box[1] + self._logo_box[3]) / 2)
            image.paste(logo, (center_x - logo.width // 2, center_y - logo.height // 2), logo)
        except Exception:
            pass
        return ImageTk.PhotoImage(image)

    def _draw_graphics_on_canvas(self):
        """
        Graphics drawn with canvas items (if Pillow is not available). Returns the x position of the title.
        """
        for row in range(0, self.height_px):
            fraction = row / (self.height_px - 1)
            color = '#%02x%02x%02x' % tuple(int(a + (b - a) * fraction) for a, b in zip(TOP_COLOR, BOTTOM_COLOR))
            self.canvas.create_line(0, row, self.width_px, row, fill=color, tags='background')
        title_x = 40
        try:
            logo = self.tk.PhotoImage(file=self._icon_path)
            factor = max(1, int(round(logo.width() / self.px(72))))
            self._logo = logo.subsample(factor, factor)
            self.canvas.create_oval(*[self.px(v) for v in self._logo_box], fill='white', outline='')
            self.canvas.create_image(self.px((self._logo_box[0] + self._logo_box[2]) / 2), self.px((self._logo_box[1] + self._logo_box[3]) / 2),
                                     image=self._logo, anchor='center')
            title_x = self._logo_box[2] + 18
        except Exception:
            pass
        left, top, right, bottom = self.plot
        self.canvas.create_line(self.px(left), self.px(bottom), self.px(right), self.px(bottom), fill=MUTED_COLOR)
        self.canvas.create_line(self.px(left), self.px(top), self.px(left), self.px(bottom), fill=MUTED_COLOR)
        points = [self.px(v) for i in range(len(self.angles)) for v in self._point(i)]
        self.canvas.create_line(*points, fill=CURVE_COLOR, width=self.px(3), capstyle='round', joinstyle='round')
        x, y = self._point(self.SPR_index)
        self.canvas.create_line(self.px(x), self.px(top), self.px(x), self.px(bottom), fill=ACCENT_COLOR, dash=(4, 3), width=self.px(2))
        self.canvas.create_oval(self.px(x - 5), self.px(y - 5), self.px(x + 5), self.px(y + 5), fill=ACCENT_COLOR, outline='')
        if self.TIR_index > 0:
            x_TIR, _ = self._point(self.TIR_index)
            self.canvas.create_line(self.px(x_TIR), self.px(top), self.px(x_TIR), self.px(bottom), fill=MUTED_COLOR, dash=(2, 4))
        return title_x

    # Showing and closing ---------------------------------------------------------------------------------------------
    def animate(self, duration=1.0):
        """
        Shows the splash screen for a moment before loading continues (blocking).
        """
        if self.closed:
            return
        try:
            end = time.time() + duration
            while time.time() < end and not self.closed:
                self.root.update()
                time.sleep(0.03)
        except Exception:
            self.close(fade=False)

    def animate_in_background(self, duration=1.0):
        """
        Kept for compatibility (the SPR curve is drawn at once).
        """
        return

    def status(self, text):
        if self.closed:
            return
        try:
            self.canvas.itemconfigure(self.status_item, text=text if len(text) <= 84 else text[:81] + '...')
            self.root.update()
        except Exception:
            pass

    def close(self, fade=True):
        if self.closed:
            return
        self.closed = True
        try:
            if getattr(self, '_animation_timer', None):
                self.root.after_cancel(self._animation_timer)
            if fade:
                for alpha in np.linspace(1, 0, 12):
                    self.root.attributes('-alpha', float(alpha))
                    self.root.update()
                    time.sleep(0.02)
            self.root.destroy()
        except Exception:
            pass


class StartMenu(SplashScreen):
    """
    Start menu with three choices:
      1) Open a measurement: a converted .csv file, or a .spr2 file that is converted first (then the channel is chosen)
      2) Reopen a previous session
      3) Convert .spr2 files to .csv without starting the rest of SPRpy (the menu stays open)
    Closing the menu (x or Escape) quits SPRpy.

    run() returns {'action': 'measurement' or 'session', 'path': ...}, or None if the menu was closed. After a choice,
    the window stays open and is used as the loading screen (status() and close()).
    """

    def __init__(self, version, sensor_values=None, data_folder=None, session_folder=None):
        super().__init__(version, sensor_values, height=MENU_HEIGHT, topmost=True)  # Always on top (its file dialogs too)
        self.data_folder = data_folder
        self.session_folder = session_folder
        self.result = None
        self.polycoeffs = None  # X-calibration of the instrument, reused for several conversions

        # Close button and keyboard shortcut
        close_item = self.text(WIDTH - 18, 14, '✕', size=18, color=MUTED_COLOR, tags='close')
        self.canvas.tag_bind(close_item, '<Button-1>', lambda event: self._quit())
        self.canvas.tag_bind(close_item, '<Enter>', lambda event: self.canvas.itemconfigure(close_item, fill=TEXT_COLOR))
        self.canvas.tag_bind(close_item, '<Leave>', lambda event: self.canvas.itemconfigure(close_item, fill=MUTED_COLOR))
        self.root.bind('<Escape>', lambda event: self._quit())

        # The window has no title bar, so it can be moved by dragging the background
        self.canvas.tag_bind('background', '<ButtonPress-1>', self._start_move)
        self.canvas.tag_bind('background', '<B1-Motion>', self._move)

        self._show_main_buttons()
        self.status('Choose how to start, or close this window (Esc) to quit.')

    # Window handling -------------------------------------------------------------------------------------------------
    def _start_move(self, event):
        self._move_offset = (event.x, event.y)

    def _move(self, event):
        x = self.root.winfo_pointerx() - self._move_offset[0]
        y = self.root.winfo_pointery() - self._move_offset[1]
        self.root.geometry('+{x}+{y}'.format(x=x, y=y))

    def _quit(self):
        self.result = None
        self.close(fade=False)

    def _finish(self, action, path):
        self.result = {'action': action, 'path': path.replace('\\', '/')}
        self._clear_buttons()
        self.root.attributes('-topmost', True)
        self.status('Loading {name}...'.format(name=os.path.basename(path)))
        self.root.quit()

    def _clear_buttons(self):
        for widget in getattr(self, '_buttons', []):
            widget.destroy()
        self._buttons = []
        self.canvas.delete('menu')

    def _button(self, text, command, x, y, width, height, primary=False):
        button = self.tk.Button(self.root, text=text, command=command, cursor='hand2', relief='flat', bd=0,
                                bg=ACCENT_COLOR if primary else BUTTON_COLOR, fg=BUTTON_TEXT_COLOR,
                                activebackground=BUTTON_ACTIVE_COLOR, activeforeground=BUTTON_TEXT_COLOR,
                                font=self.font(14, 'bold'), justify='center', wraplength=self.px(width - 16))
        self.canvas.create_window(self.px(x), self.px(y), window=button, width=self.px(width), height=self.px(height), anchor='nw', tags='menu')
        self._buttons.append(button)
        return button

    def _show_main_buttons(self):
        self._clear_buttons()
        width, height, gap = 180, 64, 20
        x0 = (WIDTH - 3 * width - 2 * gap) / 2
        y = 330
        self._button('Open measurement\n(.csv or .spr2)', self._open_measurement, x0, y, width, height, primary=True)
        self._button('Reopen session', self._open_session, x0 + width + gap, y, width, height)
        self._button('Convert .spr2 files\nonly', self._convert_only, x0 + 2 * (width + gap), y, width, height)

    def _show_channel_buttons(self, csv_files):
        """
        Buttons to choose which channel of a converted .spr2 file to open.
        """
        self._clear_buttons()
        self.text(WIDTH / 2, 332, 'Which channel do you want to open?', size=15, weight='bold', tags='menu')
        count = len(csv_files) + 1
        gap = 10
        width = min(120, (WIDTH - 60 - gap * (count - 1)) / count)
        x0 = (WIDTH - count * width - (count - 1) * gap) / 2
        for index, csv_file in enumerate(csv_files):
            label = os.path.basename(csv_file)[:-4].split('-')[-1].replace('_', '\n').replace('nm', ' nm')
            self._button(label, lambda path=csv_file: self._finish('measurement', path), x0 + index * (width + gap), 352, width, 54, primary=True)
        self._button('Back', self._back, x0 + len(csv_files) * (width + gap), 352, width, 54)

    def _back(self):
        self._show_main_buttons()
        self.status('Choose how to start, or close this window (Esc) to quit.')

    # File dialogs ----------------------------------------------------------------------------------------------------
    def _ask_file(self, title, file_types, folder):
        from tkinter.filedialog import askopenfilename
        return askopenfilename(parent=self.root, title=title, filetypes=file_types, initialdir=folder)

    def _ask_files(self, title, file_types, folder):
        from tkinter.filedialog import askopenfilenames
        return askopenfilenames(parent=self.root, title=title, filetypes=file_types, initialdir=folder)

    def _set_folder(self, path):
        if path:
            self.data_folder = os.path.dirname(path)

    # Actions ---------------------------------------------------------------------------------------------------------
    def _open_measurement(self):
        path = self._ask_file('Select a measurement file (.csv, or .spr2 to convert it first)',
                              [('SPR measurements', '*.csv *.spr2'), ('Converted measurements', '*.csv'), ('BioNavis measurements', '*.spr2')],
                              self.data_folder)
        if not path:
            self.status('No file selected. Choose how to start, or close this window (Esc) to quit.')
            return
        self._set_folder(path)
        if path.lower().endswith('.csv'):
            self._finish('measurement', path)
            return

        csv_files = self._convert([path])
        if csv_files:
            if len(csv_files) == 1:
                self._finish('measurement', csv_files[0])
            else:
                self._show_channel_buttons(csv_files)
                self.status('Converted {name} into {n} files (one per channel).'.format(name=os.path.basename(path), n=len(csv_files)))

    def _open_session(self):
        path = self._ask_file('Select a previous session file', [('SPRpy session files', '*.pickle')], self.session_folder)
        if not path:
            self.status('No session selected. Choose how to start, or close this window (Esc) to quit.')
            return
        self._finish('session', path)

    def _convert_only(self):
        paths = self._ask_files('Select .spr2 files to convert', [('BioNavis measurements', '*.spr2')], self.data_folder)
        if not paths:
            self.status('No files selected. Choose how to start, or close this window (Esc) to quit.')
            return
        self._set_folder(paths[0])
        csv_files = self._convert(paths)
        if csv_files:
            self.status('Done! Converted {n} file(s) into {m} .csv files in {folder}'.format(
                n=len(paths), m=len(csv_files), folder=os.path.dirname(csv_files[0])))

    def _convert(self, spr2_files):
        """
        Converts .spr2 files (the menu shows the progress). Returns the created .csv files, or None if it failed.
        """
        for button in self._buttons:
            button.configure(state='disabled')
        self.root.configure(cursor='watch')
        csv_files = []
        try:
            import SPRpy_spr2_to_csv
            for spr2_file in spr2_files:
                started = time.time()
                self.polycoeffs = SPRpy_spr2_to_csv.convert_spr2_file(spr2_file, self.polycoeffs, progress=self.status)
                head = os.path.join(os.path.dirname(spr2_file), os.path.basename(spr2_file)[:-5])
                created = [f for f in glob.glob(glob.escape(head) + '-L*_*nm.csv') if os.path.getmtime(f) >= started - 1]
                csv_files.extend(sorted(created, key=lambda f: int(os.path.basename(f)[len(os.path.basename(head)) + 2:].split('_')[0])))
            return csv_files
        except Exception as error:
            traceback.print_exc()
            self.status('Conversion failed: {error}'.format(error=error))
            return None
        finally:
            if not self.closed:
                self.root.configure(cursor='')
                for button in self._buttons:
                    try:
                        button.configure(state='normal')
                    except Exception:
                        pass

    def run(self):
        """
        Shows the menu until a measurement or session is chosen (returns the choice) or the menu is closed (None).
        """
        self.root.lift()
        self.root.focus_force()
        self.animate_in_background()
        self.root.mainloop()
        return None if self.closed else self.result


def show_splash(version, sensor_values=None, animate=True):
    """
    Creates and shows the splash screen. Returns None if it can not be shown (e.g. no display).
    """
    try:
        splash = SplashScreen(version, sensor_values)
        if animate:
            splash.animate()
        return splash
    except Exception:
        return None


def open_start_menu(version, sensor_values=None, data_folder=None, session_folder=None):
    """
    Creates the start menu. Returns None if it can not be shown (e.g. no display), so that the console prompts are used.
    """
    try:
        return StartMenu(version, sensor_values, data_folder=data_folder, session_folder=session_folder)
    except Exception:
        traceback.print_exc()
        return None
