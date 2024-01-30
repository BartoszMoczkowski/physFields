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



class CustomDropdown(DropDown):
    def __init__(self, **kwargs):
        
        super().__init__(**kwargs)
        self.select('item1')


class DisplayCanvas(BoxLayout):
    def __init__(self, **kwargs):


        self.elements = 1 # 1 - rectangle , 2 - traingle
        self.state = 1 # 1 - add material 2 - remove material 3 - add heat sourcess 4 remove heat sources
        self.shape = (100,100)
        self.heatmap = np.zeros(self.shape)
        self.is_drawing = False 
        Window.bind(on_motion = self.on_motion)
        super().__init__(**kwargs)


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

        if self.state == 1:
            self.heatmap[x_0:x_n+1,y_0:y_n+1] = 1
        elif self.state == 2:
            self.heatmap[x_0:x_n+1,y_0:y_n+1] = 0
        self.update_display()


    def update_display(self):
        Rect = self
        rect_w = int(Rect.width)
        rect_h = int(Rect.height)
        texture_shape = self.shape
        texture = Texture.create(size=texture_shape)

        buf = (np.repeat(self.heatmap.flatten(),4)*255).astype(np.uint8)
        
        texture.blit_buffer(buf, colorfmt='rgba', bufferfmt='ubyte')
        Rect.canvas.clear()
        with Rect.canvas:
            Rectangle(texture=texture, pos=Rect.pos, size=(rect_w, rect_h))

    def button_set_state(self,state):
        self.state = state
class Menu(BoxLayout):

    def __init__(self, **kwargs):
        
        Window.bind(on_resize=self.on_resize)
        super().__init__(**kwargs)

    def on_resize(self,*args):
        self.onApply()

    def onApply(self):
        Rect = self.ids.Display
        Rect.update_display()


       
class MyApp(App):

    def __init__(self, **kwargs):

        super().__init__(**kwargs)
    def build(self):

        return Menu()
    

if __name__ == '__main__':
    MyApp().run()