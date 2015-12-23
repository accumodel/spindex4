#-------------------------------------------------------------------------------
# Name:        spindex4
# Purpose:     a simple quadtree spatial index (not cythonized version)
#
# Author:      Mark Wilson
# Date:        2015-12-23
# Notes:
#
# Licence:     MIT
#-------------------------------------------------------------------------------

# ******************************************************************************
# Quad Tree implementation for points (SpIndex4p) and for envelopes (SpIndex4e)
#
# It anticipates that the user will seed the index with the approximate bounds
# of all of the features but it also has an expand method if an inserted item
# falls outside of the current bounds.
#
# Currently negative coordinates are not supported
# ******************************************************************************

import copy

class SpIndex4(object):
    '''quad index superclass with an anticipated outer boundary
    assumes positive coordinates'''
    def __cinit__(self, tuple bnds):
        assert (bnds[2] > bnds[0] and bnds[3] > bnds[1]), 'The boundary object is not a proper envelope'
        cdef tuple bounds
        for item in bnds:
            assert(item >= 0), 'The bounds are in the negative domain'
        self.bounds = (float(bnds[0]),  # use float in case the envelope coordinates are passed in as strings
                       float(bnds[1]),
                       float(bnds[2]),
                       float(bnds[3]))
    def _intersect(self,tuple envelope1,tuple envelope2):
        '''returns True if the two envelopes overlap'''
        return (envelope1[0] <= envelope2[2] and
                envelope1[1] <= envelope2[3] and
                envelope1[2] >= envelope2[0] and
                envelope1[3] >= envelope2[1])
    def _point_in_envelope(self,tuple point,tuple envelope):
        return (point[0] <= envelope[2] and
                point[1] <= envelope[3] and
                point[0] >= envelope[0] and
                point[1] >= envelope[1])
    def _can_fit_in_subcell(self,tuple obj_env,tuple cell_env):
        x_mid = (cell_env[2] + cell_env[0]) / 2.0
        y_mid = (cell_env[3] + cell_env[1]) / 2.0
        return (self._contains(obj_env,(cell_env[0],cell_env[1],x_mid,y_mid)) or
                self._contains(obj_env,(x_mid,cell_env[1],cell_env[2],y_mid)) or
                self._contains(obj_env,(cell_env[0],y_mid,x_mid,cell_env[3])) or
                self._contains(obj_env,(x_mid,y_mid,cell_env[2],cell_env[3])))


class SpIndex4p(SpIndex4):
    '''Uses SpIndex4 superclass and implements a point specific index
    the point index stores the point instead of an envelope'''
    def __init__(self,bounds = (0,0,10000,10000)):
        '''initializer that sets the initial bounds then breaks it into four
        quadrants/levels via the _divide function'''
        SpIndex4.__cinit__(self,bounds)
        self._index = self._divide(self.bounds)
    def _divide(self, tuple envelope):
        '''divides one level into four sub levels'''
        cdef double minx = envelope[0]
        cdef double miny = envelope[1]
        cdef double maxx = envelope[2]
        cdef double maxy = envelope[3]
        cdef double x_sp = (maxx - minx) / 2.0
        cdef double y_sp = (maxy - miny) / 2.0
        return [[(minx,miny,minx + x_sp, miny + y_sp),[]],
                [(minx + x_sp,miny,maxx, miny + y_sp),[]],
                [(minx,miny + y_sp,minx + x_sp, maxy),[]],
                [(minx + x_sp,miny + y_sp,maxx, maxy),[]]]
    def _expand(self, tuple pt):
        '''if the user attempts an insertion outside the current
        bounds, the grid will be automatically expanded'''
        cdef double x_midpt = (self.bounds[2] + self.bounds[0])/2.0
        cdef double y_midpt = (self.bounds[3] + self.bounds[1])/2.0
        cdef double old_xmin = self.bounds[0]
        cdef double old_ymin = self.bounds[1]
        cdef double old_xmax = self.bounds[2]
        cdef double old_ymax = self.bounds[3]
        cdef double old_width = old_xmax - old_xmin
        cdef double old_height = old_ymax - old_ymin
        cdef tuple old_bounds = self.bounds
        cdef double new_xmin
        cdef double new_ymin
        if pt[0] < x_midpt:     # lower x bounds
           if pt[1] < y_midpt:  # lower x bounds lower y bounds expansion
              new_xmin = old_xmin - old_width
              new_ymin = old_ymin - old_height
              self.bounds = (new_xmin,
                             new_ymin,
                             old_xmax,
                             old_ymax)
              self._index = [[(new_xmin,new_ymin,old_xmin,old_ymin),[]],
                             [(old_xmin,new_ymin,old_xmax,old_ymin),[]],
                             [(new_xmin,old_ymin,old_xmin,old_ymax),[]],
                             [old_bounds,self._index]]  #[0]]
           else:  # lower x bounds upper y bounds expansion
              new_xmin = old_xmin - old_width
              new_ymax = old_ymax + old_height
              self.bounds = (new_xmin,
                             old_ymin,
                             old_xmax,
                             new_ymax)
              self._index = [[(new_xmin,old_ymin,old_xmin,old_ymax),[]],
                             [old_bounds,self._index],  #[0],
                             [(new_xmin,old_ymax,old_xmin,new_ymax),[]],
                             [(old_xmin,old_ymax,old_xmax,new_ymax),[]]]
        else:      # upper x bounds
           if pt[1] < y_midpt:  # upper x bounds lower y bounds expansion
              new_ymin = old_ymin - old_height
              new_xmax = old_xmax + old_width
              self.bounds = (old_xmin,
                             new_ymin,
                             new_xmax,
                             old_ymax)
              self._index = [[(old_xmin,new_ymin,old_xmax,old_ymin),[]],
                             [(old_xmax,new_ymin,new_xmax,old_ymin),[]],
                             [old_bounds,self._index],  #[0],
                             [(old_xmax,old_ymin,new_xmax,old_ymax),[]]]
           else:  # upper x bounds upper y bounds expansion
              new_xmax = old_xmax + old_width
              new_ymax = old_ymax + old_height
              self.bounds = (old_xmin,
                             old_ymin,
                             new_xmax,
                             new_ymax)
              self._index = [[old_bounds,self._index],  #[0],
                             [(old_xmax,old_ymin,new_xmax,old_ymax),[]],
                             [(old_xmin,old_ymax,old_xmax,new_ymax),[]],
                             [(old_xmax,old_ymax,new_xmax,new_ymax),[]]]
        if not self._contains(pt,self.bounds):
           self._expand(pt)
    def insert(self, tuple object_coord, object_key, leaf = None):
        '''inserts an X,Y coordinate into the index.  if there is already a point
        in the leaf, then it subdivides (unless the coordinates are exactly the same)'''
        if not self._contains(object_coord,self.bounds): # The object you are trying to insert is outside the bounds of the spatial index
           self._expand(object_coord)
        if not leaf:  # set current leaf to main index if this is not a recursive call
           leaf = self._index
        for cell in leaf:
            if self._contains(object_coord,cell[0]):
               if len(cell[1]) == 0:  # empty cell
                  if type(object_key) != set:
                     object_key = set([object_key])
                  cell[1] = [object_key,object_coord]  # insert the object
               elif len(cell[1]) == 2:  # this is an item
                  if object_coord == cell[1][1]:  # object at the same envelope
                     cell[1][0].add(object_key)  # so just append to the list
                  else:
                     tempID = copy.deepcopy(cell[1][0])  # copy the item and insert it again at a lower level
                     temp_coord = cell[1][1]
                     cell[1] = self._divide(cell[0])
                     # recursive calls
                     self.insert(temp_coord,tempID,cell[1])  # push the item down to the newly divided cell
                     self.insert(object_coord,object_key,cell[1])  # try to insert the original item
               else:
                  self.insert(object_coord,object_key,cell[1])  # try to insert the original item
    def remove(self, object_key, leaf = None):  # TODO.  Verify functionality.  MAY BE TOO COMPLEX
        '''searches for and removes the object from the index'''
        if not leaf:
           leaf = self._index
        for cell in leaf:
            if cell[1]:  # not empty
              if len(cell[1]) > 2:  # recursive
                 for item in self.remove(object_key, cell[1]):
                     yield item[0],item[1]
              else:
                 if object_key in cell[1][0]:
                     if len(cell[1][0]) == 1:
                         pass
    def intersect(self, envelope, leaf = None):
        '''returns all IDs with coordinates that fall inside the envelope that is passed in'''
        if not leaf:
           leaf = self._index
        for cell in leaf:
            cell_env = cell[0]
            cell_contents = cell[1]
            if len(cell_contents) > 0 and self._intersect(envelope,cell_env):
               if len(cell_contents) == 2:  # an item
                  itemIDs = cell_contents[0]
                  item_coord = cell_contents[1]
                  if self._contains(item_coord,envelope):
                     for item in itemIDs:
                         yield item
               else:
                  # recursive call
                  for item in self.intersect(envelope,cell_contents):
                      yield item
    def _contains(self,tuple point, tuple envelope):
        '''envelope COMPLETELY contains point'''
        return (envelope[0] <= point[0] and
                envelope[1] <= point[1] and
                envelope[2] >= point[0] and
                envelope[3] >= point[1])
    def walk(self, level = 0, leaf = None):
        '''mainly a troubleshooting function to return all of the leafs'''
        if not leaf:
           leaf = self._index
        for cell in leaf:
            if cell[1]:  # not empty
              if len(cell[1]) > 2:  # recursive
                 for item in self.walk(level + 1, cell[1]):
                     yield item[0],item[1]
              else:
                 yield cell[1],level
    def rectangles(self,leaf = None):
        '''returns all of the grid rectangles so that they can be plotted'''
        if not leaf:
           leaf = self._index
        for cell in leaf:
            yield cell[0]
            if len(cell[1]) == 4:
               for item in self.rectangles(cell[1]):
                   yield item
    def print_rectangles(self,level = 0, leaf = None):
        '''mainly troubleshooting function'''
        if not leaf:
           leaf = self._index
        for cell in leaf:
            if len(cell[1]) == 4:
               print '->' * level, cell[0], ' see subleaf'
               self.print_rectangles(level + 1, cell[1])
            else:
                 print '->' * level, cell
    def visualize(self):
        '''troubleshooting function that plots stuff on screen (via Tk) so it
        can be visualized'''
        # set scale so that maxx or maxy is 1000
        from Tkinter import Tk, Canvas, mainloop, NW
        master = Tk()
        w = 1000
        h = 1000
        wnd = Canvas(master, width=w, height=h)
        wnd.xview_moveto(0)
        wnd.yview_moveto(0)
        wnd.pack()
        scale = min(w/self.bounds[2], h/self.bounds[3])
        for item in self.rectangles():
            wnd.create_rectangle([val * scale for val in item])
        for item in self.walk():
            pt = (item[0][1][0]*scale-1,
                  item[0][1][1]*scale-1,
                  item[0][1][0]*scale+1,
                  item[0][1][1]*scale+1)
            wnd.create_rectangle(pt, fill='red')
            wnd.create_text((item[0][1][0]*scale,item[0][1][1]*scale), text=','.join(item[0][0]), anchor=NW)
        mainloop()
        del w, master

class SpIndex4e(SpIndex4):
    '''Uses SpIndex4 superclass and implements an envelope (line, polygon) specific index'''
    def __cinit__(self,tuple bounds = (0,0,10000,10000)):
        SpIndex4.__cinit__(self,bounds)
        self._index = self._divide(self.bounds)
    def _divide(self, tuple envelope):
        '''divides one level into four sub levels'''
        cdef double minx = envelope[0]
        cdef double miny = envelope[1]
        cdef double maxx = envelope[2]
        cdef double maxy = envelope[3]
        cdef double x_sp = (maxx - minx) / 2.0
        cdef double y_sp = (maxy - miny) / 2.0
        return [[(minx,miny,minx + x_sp, miny + y_sp),{}],
                [(minx + x_sp,miny,maxx, miny + y_sp),{}],
                [(minx,miny + y_sp,minx + x_sp, maxy),{}],
                [(minx + x_sp,miny + y_sp,maxx, maxy),{}]]
    def _expand(self, tuple env):
        '''if the user attempts an insertion outside the current
        bounds, the grid will be automatically expanded'''
        cdef double new_xmin
        cdef double new_ymin
        cdef double env_cent_x = (env[2] + env[0]) / 2.0
        cdef double env_cent_y = (env[3] + env[1]) / 2.0
        pt = env_cent_x, env_cent_y
        cdef double x_midpt = (self.bounds[2] + self.bounds[0]) / 2.0
        cdef double y_midpt = (self.bounds[3] + self.bounds[1]) / 2.0
        cdef double old_xmin = self.bounds[0]
        cdef double old_ymin = self.bounds[1]
        cdef double old_xmax = self.bounds[2]
        cdef double old_ymax = self.bounds[3]
        cdef double old_width = old_xmax - old_xmin
        cdef double old_height = old_ymax - old_ymin
        cdef tuple old_bounds = self.bounds
        if pt[0] < x_midpt:     # lower x bounds
           if pt[1] < y_midpt:  # lower x bounds lower y bounds expansion
              new_xmin = old_xmin - old_width
              new_ymin = old_ymin - old_height
              self.bounds = (new_xmin,
                             new_ymin,
                             old_xmax,
                             old_ymax)
              self._index = [[(new_xmin,new_ymin,old_xmin,old_ymin),{}],
                             [(old_xmin,new_ymin,old_xmax,old_ymin),{}],
                             [(new_xmin,old_ymin,old_xmin,old_ymax),{}],
                             [old_bounds,self._index]]  #[0]]
           else:  # lower x bounds upper y bounds expansion
              new_xmin = old_xmin - old_width
              new_ymax = old_ymax + old_height
              self.bounds = (new_xmin,
                             old_ymin,
                             old_xmax,
                             new_ymax)
              self._index = [[(new_xmin,old_ymin,old_xmin,old_ymax),{}],
                             [old_bounds,self._index],  #[0],
                             [(new_xmin,old_ymax,old_xmin,new_ymax),{}],
                             [(old_xmin,old_ymax,old_xmax,new_ymax),{}]]
        else:      # upper x bounds
           if pt[1] < y_midpt:  # upper x bounds lower y bounds expansion
              new_ymin = old_ymin - old_height
              new_xmax = old_xmax + old_width
              self.bounds = (old_xmin,
                             new_ymin,
                             new_xmax,
                             old_ymax)
              self._index = [[(old_xmin,new_ymin,old_xmax,old_ymin),{}],
                             [(old_xmax,new_ymin,new_xmax,old_ymin),{}],
                             [old_bounds,self._index],  #[0],
                             [(old_xmax,old_ymin,new_xmax,old_ymax),{}]]
           else:  # upper x bounds upper y bounds expansion
              new_xmax = old_xmax + old_width
              new_ymax = old_ymax + old_height
              self.bounds = (old_xmin,
                             old_ymin,
                             new_xmax,
                             new_ymax)
              self._index = [[old_bounds,self._index],  #[0],
                             [(old_xmax,old_ymin,new_xmax,old_ymax),{}],
                             [(old_xmin,old_ymax,old_xmax,new_ymax),{}],
                             [(old_xmax,old_ymax,new_xmax,new_ymax),{}]]
        if not self._contains(env,self.bounds):
           self._expand(env)
    def insert(self, tuple object_envelope, object_key, leaf = None):
        '''adds the object ID to any leaf that it intersects.  it only subdivides
        the leaf if the object envelope is small enough to fit inside a quadrant if
        the leaf/cell was subdivided'''
        if not self._contains(object_envelope,self.bounds):
           self._expand(object_envelope)
        if not leaf:  # if this is not a recursive call then set the leaf to the whole index
           leaf = self._index
        for cell in leaf:
            if self._intersect(object_envelope,cell[0]):
               if type(cell[1]) == list:  # there is a subleaf
                  self.insert(object_envelope,object_key, cell[1])
               else:
                   # if the object can fit completely into a subcell then divide
                   if self._can_fit_in_subcell(object_envelope,cell[0]):
                      temp_dict = cell[1]
                      cell[1] = self._divide(cell[0])
                      # reinsert each existing item
                      for key in temp_dict:
                         self.insert(temp_dict[key],key,cell[1])
                      # now insert the new item
                      self.insert(object_envelope,object_key,cell[1])
                   # otherwise, just insert it
                   else:
                      cell[1][object_key] = object_envelope
    def intersect(self, tuple envelope, leaf = None):
        '''returns all object IDs if their envelope intersects the passed in envelope'''
        if not leaf:
           leaf = self._index
        for cell in leaf:
            cell_env = cell[0]
            cell_contents = cell[1]
            if len(cell_contents) > 0 and self._intersect(envelope,cell_env):
               if type(cell_contents) == dict:  # an item
                  for key in cell_contents:
                      if self._intersect(envelope,cell_contents[key]):
                         yield key
               else:
                  for item in self.intersect(envelope,cell_contents):
                      yield item
    def _contains(self,tuple envelope1,tuple envelope2):
        '''envelope2 COMPLETELY contains envelope1'''
        return (envelope1[0] >= envelope2[0] and
                envelope1[1] >= envelope2[1] and
                envelope1[2] <= envelope2[2] and
                envelope1[3] <= envelope2[3])
    def walk(self, level = 0, leaf = None):
        '''mainly a troubleshooting function'''
        if not leaf:
           leaf = self._index
        for cell in leaf:
            if len(cell[1]) > 0:  # not empty
              if type(cell[1]) == list:  # recursive
                 for item in self.walk(level + 1, cell[1]):
                     yield item[0],item[1]
              else:
                 yield cell[1],level
    def rectangles(self,leaf = None):
        '''returns all of the grid rectangles so that they can be plotted'''
        if not leaf:
           leaf = self._index
        for cell in leaf:
            yield cell[0]
            if type(cell[1]) == list:
               for item in self.rectangles(cell[1]):
                   yield item
    def print_rectangles(self,level = 0, leaf = None):
        '''mainly a troubleshooting function'''
        if not leaf:
           leaf = self._index
        for cell in leaf:
            if type(cell[1]) == list:
               print '->' * level, cell[0], ' see subleaf'
               self.print_rectangles(level + 1, cell[1])
            else:
                 print '->' * level, cell
    def visualize(self):
        '''mainly a troubleshooting function'''
        # set scale so that maxx or maxy is 1000
        from Tkinter import  Tk, Canvas, mainloop, NW
        master = Tk()
        w = 1000
        h = 1000
        wnd = Canvas(master, width=w, height=h)
        wnd.xview_moveto(0)
        wnd.yview_moveto(0)
        wnd.pack()
        scale = min(w/self.bounds[2], h/self.bounds[3])
        for item in self.rectangles():
            wnd.create_rectangle([val * scale for val in item])
        for item in self.walk():
            for key in item[0]:
                rect = (item[0][key][0]*scale,
                        item[0][key][1]*scale,
                        item[0][key][2]*scale,
                        item[0][key][3]*scale)
                wnd.create_rectangle(rect)
                wnd.create_text((item[0][key][0]*scale,item[0][key][1]*scale), text=key, anchor=NW)
        mainloop()
        del w, master


def main():
    pass

if __name__ == '__main__':
    main()


