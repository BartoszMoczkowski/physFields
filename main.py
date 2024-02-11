import kivy
kivy.require('2.1.0')

#biblioteki to kivy i numpy

from kivy.app import App
from kivy.uix.image import Image
from kivy.graphics.texture import Texture
from kivy.graphics import Rectangle,Color
from kivy.uix.boxlayout import BoxLayout
import numpy as np
from kivy.uix.widget import Widget
from kivy.core.window import Window
from kivy.uix.dropdown import DropDown
from kivy.uix.textinput import TextInput
import matplotlib.colors as cl
from FEM_module import FemSolver

class IntInput(TextInput):

    def insert_text(self, substring, from_undo=False):

        try:
            int(substring)
        except:
            substring = ""
        return super().insert_text(substring, from_undo)


class FloatInput(TextInput):
    def insert_text(self, substring, from_undo=False):
        try:
         
            float(self.text + substring)
        except:
            if substring != "-":
                substring = ""
            else:
                substring = "-"
        return super().insert_text(substring, from_undo)
class CustomDropdown(DropDown):
    def __init__(self, **kwargs):
        
        super().__init__(**kwargs)
        self.select('item1')


class DisplayCanvas(BoxLayout):
    def __init__(self, **kwargs):

        self.fem = FemSolver()
        self.fem.n_element_x = 20
        self.fem.n_element_y = 10

        self.fem.dirchlet_conditions[0, :] = 10
        self.fem.dirchlet_conditions[-1,:] = 100

        self.elements = 1 # 1 - rectangle , 2 - traingle
        self.editing_state = 1 # 1 - add material 2 - remove material 
        self.layer_select = 0 #0 - nodes 1 - material 2 - heat source
        self.layer_selected = False # 

        self.shape = (self.fem.n_element_y,self.fem.n_element_x)
        self.heatmap = np.ones(self.shape)
        
        
        self.is_drawing = False 
        Window.bind(on_motion = self.on_motion)
        super().__init__(**kwargs)


    def get_root_app(self):
        if App.get_running_app() != None:
            return App.get_running_app()            
        
    def get_ids(self):
        return self.get_root_app().root.ids

    def on_touch_down(self, touch):
        if self.collide_point(*touch.pos):
            self.is_drawing = True
            self.menu = App.get_running_app().root
            self.last_click_pos = touch.pos
            offset = np.array(touch.pos)-np.array(self.pos)
            offset_w = offset[0]/self.width
            offset_h = offset[1]/self.height
            self.pixel_offset_w_start = int(offset_w * self.shape[0])
            self.pixel_offset_h_start = int(offset_h * self.shape[1]) 
            

            touch.grab(self)
        return super().on_touch_down(touch)
    
    def on_touch_up(self, touch):
        if self.collide_point(*touch.pos) and self.is_drawing == True:

            self.is_drawing = False
            offset = np.array(touch.pos)-np.array(self.pos)
            offset_w = offset[0]/self.width
            offset_h = offset[1]/self.height

            
            self.pixel_offset_w_end = int(offset_w * self.shape[0])
            self.pixel_offset_h_end = int(offset_h * self.shape[1]) 

            self.change_array() 

            touch.ungrab(self)
        return super().on_touch_up(touch)
    
    def on_motion(self, etype, me,x):
        if self.is_drawing and not self.collide_point(*x.pos):
            self.update_display()
            return
        if self.is_drawing:
            if self.last_click_pos[0] < x.pos[0]:
                x_0 = self.last_click_pos[0]
                x_n = x.pos[0]
            else:
                x_0 = x.pos[0]  
                x_n = self.last_click_pos[0]

            if self.last_click_pos[1] < x.pos[1]:
                y_0 = self.last_click_pos[1]
                y_n = x.pos[1]
            else:
                y_0 = x.pos[1]  
                y_n = self.last_click_pos[1]
            
            self.update_display()
            with self.canvas:
                Color(1,1,1,1.)

                Rectangle(pos = (x_0,y_0), size = (abs(self.last_click_pos[0] - x.pos[0]),abs(self.last_click_pos[1] - x.pos[1])))
                Color(0,0,0,1.)
                Rectangle(pos = (x_0+2,y_0+2), size = (abs(self.last_click_pos[0] - x.pos[0])-4,abs(self.last_click_pos[1] - x.pos[1])-4))

    def change_array(self):
        if self.pixel_offset_h_start > self.pixel_offset_h_end:
            x_0 = self.pixel_offset_h_end
            x_n = self.pixel_offset_h_start
        else:
            x_0 = self.pixel_offset_h_start
            x_n = self.pixel_offset_h_end

        if self.pixel_offset_w_start > self.pixel_offset_w_end:
            y_0 = self.pixel_offset_w_end
            y_n = self.pixel_offset_w_start
        else:
            y_0 = self.pixel_offset_w_start
            y_n = self.pixel_offset_w_end



        input_text = self.get_ids().value_input.text
        editing_value = 0
        if self.editing_state == 1 and input_text != "" :
            try:
                editing_value = float(input_text) 
                if editing_value < 0 and self.layer_select in [0,1,2]:
                    editing_value = 0
            except:
                editing_value = 0
        elif self.editing_state == 2:
            editing_value = 0

        if self.layer_select == 0:
            editing_value = 1 if editing_value != 0 else 0
            self.fem.space[x_0:x_n+1,y_0:y_n+1] = editing_value
        if self.layer_select == 1:
            self.fem.material[x_0:x_n+1,y_0:y_n+1] = editing_value
        if self.layer_select == 2:
            self.fem.dirchlet_conditions[x_0:x_n+1,y_0:y_n+1] = editing_value
        if self.layer_select == 3:
            self.fem.neumann_conditions[x_0:x_n+1,y_0:y_n+1] = editing_value

        self.update_display()


    def update_display(self):
        Rect = self
        rect_w = int(Rect.width)
        rect_h = int(Rect.height)
        texture_shape = self.shape
        texture = Texture.create(size=texture_shape)
        texture.mag_filter = "nearest"


        buf_source = self.heatmap

        if self.layer_selected:
            if self.layer_select == 0:
                buf_source = self.fem.space
            if self.layer_select == 1:
                buf_source = self.fem.material
            if self.layer_select == 2:
                buf_source = self.fem.dirchlet_conditions
            if self.layer_select == 3:
                buf_source = self.fem.neumann_conditions

        if buf_source.min() < 0:
            buf_source = buf_source + abs(buf_source.min()) + 0.01
        if buf_source.max() != 0:
            buf_source = buf_source/buf_source.max()
        h = (1 - buf_source.flatten()) * (240/360)
        s = np.ones_like(h)
        v = np.ones_like(h) - 0.3
        buf = (cl.hsv_to_rgb(np.column_stack([h,s,v])).flatten()*255).astype(np.uint8)
        
        texture.blit_buffer(buf, colorfmt='rgb', bufferfmt='ubyte')
        Rect.canvas.clear()
        with Rect.canvas:
            Rectangle(texture=texture, pos=Rect.pos, size=(rect_w, rect_h))

    def solve(self):
        self.fem.n_element_x = int(self.get_ids().width_input.text)
        self.fem.n_element_y = int(self.get_ids().height_input.text)
        self.shape = (self.fem.n_element_x,self.fem.n_element_y)
        self.heatmap = np.ones(self.shape)
        self.heatmap = self.fem.solve()
        self.heatmap = self.heatmap

        self.update_display()

    def button_set_state(self,state):
        self.editing_state = state

    def button_set_layer(self,layer,active):
        self.layer_selected = active
        self.layer_select = layer
        self.update_display()
class Menu(BoxLayout):

    def __init__(self, **kwargs):
        
        Window.bind(on_resize=self.on_resize)
        Window.bind(on_maximize=self.on_resize)
        Window.bind(on_restore=self.on_resize)
        super().__init__(**kwargs)

    def on_resize(self,*args):
        self.ids.Display.update_display()

    

    def onApply(self):
        Rect = self.ids.Display
        Rect.solve()

       
class MyApp(App):

    def __init__(self, **kwargs):

        super().__init__(**kwargs)
    def build(self):

        return Menu()
    



if __name__ == '__main__':
    MyApp().run()