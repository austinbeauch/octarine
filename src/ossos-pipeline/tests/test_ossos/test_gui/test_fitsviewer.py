__author__ = "David Rusk <drusk@uvic.ca>"

import unittest
import wx

from matplotlib.backend_bases import MouseEvent as MPLMouseEvent
from mock import Mock
from hamcrest import assert_that, equal_to, has_length, instance_of, none

from ossos.gui import fitsviewer
from ossos.gui.fitsviewer import (InteractionContext, MoveCircleState,
                                 CreateCircleState, AdjustColormapState)


class MPLFitsImageViewerTest(unittest.TestCase):
    def setUp(self):
        self.app = wx.App()
        self.rootframe = wx.Frame(None)

        self.viewer = fitsviewer.MPLFitsImageViewer(self.rootframe)

    def test_draw_one_circle(self):
        axes = self.viewer.axes

        assert_that(axes.patches, has_length(0))
        cx = 1
        cy = 2
        cr = 3
        self.viewer.draw_circle(cx, cy, cr)

        assert_that(axes.patches, has_length(1))
        circle = axes.patches[0]

        assert_that(circle.center, equal_to((cx, cy)))
        assert_that(circle.radius, equal_to(cr))

    def test_draw_second_circle_removes_first(self):
        axes = self.viewer.axes

        c1x = 1
        c1y = 2
        c1r = 3
        self.viewer.draw_circle(c1x, c1y, c1r)

        assert_that(axes.patches, has_length(1))

        c2x = 4
        c2y = 5
        c2r = 6
        self.viewer.draw_circle(c2x, c2y, c2r)

        assert_that(axes.patches, has_length(1))

        circle = axes.patches[0]

        assert_that(circle.center, equal_to((c2x, c2y)))
        assert_that(circle.radius, equal_to(c2r))


class InteractionTest(unittest.TestCase):
    def setUp(self):
        self.app = wx.App()
        self.rootframe = wx.Frame(None)
        self.viewer = fitsviewer.MPLFitsImageViewer(self.rootframe)
        self.viewer.figure = Mock()
        self.viewer.axes = Mock()

        self.interaction_context = InteractionContext(self.viewer)

    def _create_mouse_event(self, x, y, button, inaxes=True):
        event = Mock(spec=MPLMouseEvent)
        event.x = x
        event.xdata = x
        event.y = y
        event.ydata = y
        event.button = button

        if inaxes:
            event.inaxes = self.viewer.axes
        else:
            event.inaxes = Mock()  # a new, different axes

        return event

    def fire_press_event(self, x, y, button=InteractionContext.MOUSE_BUTTON_LEFT,
                         inaxes=True):
        self.interaction_context.on_press(
            self._create_mouse_event(x, y, button, inaxes))

    def fire_release_event(self, button=InteractionContext.MOUSE_BUTTON_LEFT):
        event = Mock(spec=MPLMouseEvent)
        event.button = button
        self.interaction_context.on_release(event)

    def fire_motion_event(self, x, y, inaxes=True):
        self.interaction_context.on_motion(
            self._create_mouse_event(x, y, inaxes))

    def test_state_click_in_circle(self):
        x = 10
        y = 10
        radius = 5

        self.viewer.draw_circle(x, y, radius)
        self.fire_press_event(x + 2, y + 2)
        assert_that(self.interaction_context.state, instance_of(MoveCircleState))

    def test_press_release(self):
        x = 10
        y = 10
        radius = 5

        self.viewer.draw_circle(x, y, radius)
        assert_that(not self.interaction_context.state.pressed)
        self.fire_press_event(x + 2, y + 2)
        assert_that(self.interaction_context.state.pressed)
        self.fire_release_event()
        assert_that(not self.interaction_context.state.pressed)

    def test_state_click_outside_circle(self):
        x = 10
        y = 10
        radius = 5

        self.viewer.draw_circle(x, y, radius)
        self.fire_press_event(x + 2, y + 2)
        assert_that(self.interaction_context.state, instance_of(MoveCircleState))
        self.fire_release_event()
        assert_that(self.interaction_context.state, instance_of(MoveCircleState))
        self.fire_press_event(x + 6, y + 6)
        assert_that(self.interaction_context.state, instance_of(CreateCircleState))

    def test_state_right_click(self):
        x = 10
        y = 10

        self.fire_press_event(x, y, button=InteractionContext.MOUSE_BUTTON_LEFT)
        assert_that(self.interaction_context.state, instance_of(CreateCircleState))
        self.fire_release_event(button=InteractionContext.MOUSE_BUTTON_LEFT)

        self.fire_press_event(x, y, button=InteractionContext.MOUSE_BUTTON_RIGHT)
        assert_that(self.interaction_context.state, instance_of(AdjustColormapState))
        self.fire_release_event(button=InteractionContext.MOUSE_BUTTON_RIGHT)

        self.fire_press_event(x, y, button=InteractionContext.MOUSE_BUTTON_LEFT)
        assert_that(self.interaction_context.state, instance_of(CreateCircleState))
        self.fire_release_event(button=InteractionContext.MOUSE_BUTTON_LEFT)

    def test_drag_circle(self):
        x0 = 10
        y0 = 10
        radius = 5

        xclick = x0 + 2
        yclick = y0 + 2
        dx = 10
        dy = 5

        self.viewer.draw_circle(x0, y0, radius)
        assert_that(self.interaction_context.get_circle().center, equal_to((x0, y0)))
        self.fire_press_event(xclick, yclick)

        self.fire_motion_event(xclick + dx, yclick + dy)
        assert_that(self.interaction_context.get_circle().center,
                    equal_to((x0 + dx, y0 + dy)))
        assert_that(self.interaction_context.get_circle().radius, equal_to(radius))

    def test_create_circle(self):
        x0 = 10
        y0 = 10
        dx = 10
        dy = 30

        assert_that(self.interaction_context.get_circle(), none())
        self.fire_press_event(x0, y0)
        self.fire_motion_event(x0 + dx, y0 + dy)
        assert_that(self.interaction_context.get_circle().center,
                    equal_to((15, 25)))
        assert_that(self.interaction_context.get_circle().radius, equal_to(15))

    def test_has_had_interaction(self):
        assert_that(not self.viewer.has_had_interaction())
        self.fire_press_event(10, 10)
        # NOTE: don't count an interaction until it has altered the viewer in
        # a visible way
        assert_that(not self.viewer.has_had_interaction())
        self.fire_motion_event(12, 12)
        assert_that(self.viewer.has_had_interaction())

    def test_motion_not_pressed(self):
        x = 10
        y = 10
        radius = 5

        self.viewer.draw_circle(x, y, radius)

        self.interaction_context.state = CreateCircleState(self.interaction_context)
        self.fire_motion_event(x + 2, y + 2)
        assert_that(self.interaction_context.get_circle().center, equal_to((x, y)))
        assert_that(self.interaction_context.get_circle().radius, equal_to(radius))

        self.interaction_context.state = MoveCircleState(self.interaction_context)
        self.fire_motion_event(x + 2, y + 2)
        assert_that(self.interaction_context.get_circle().center, equal_to((x, y)))
        assert_that(self.interaction_context.get_circle().radius, equal_to(radius))

    def test_click_no_drag_inside_circle(self):
        x = 10
        y = 10
        radius = 5

        self.viewer.draw_circle(x, y, radius)

        click_x = 12
        click_y = 13
        self.fire_press_event(click_x, click_y)
        self.fire_release_event()

        assert_that(self.interaction_context.get_circle().center,
                    equal_to((click_x, click_y)))

    def test_click_no_drag_outside_circle(self):
        x = 10
        y = 10
        radius = 5

        self.viewer.draw_circle(x, y, radius)

        click_x = 20
        click_y = 21
        self.fire_press_event(click_x, click_y)
        self.fire_release_event()

        assert_that(self.interaction_context.get_circle().center,
                    equal_to((click_x, click_y)))

    def test_xy_changed_event_on_click(self):
        handler = Mock()
        self.viewer.register_xy_changed_event_handler(handler)

        self.viewer.draw_circle(10, 10, 5)

        x_click = 20
        y_click = 30
        self.fire_press_event(x_click, y_click)
        self.fire_release_event()

        handler.assert_called_once_with(x_click, y_click)

    def test_xy_changed_event_on_drag(self):
        handler = Mock()
        self.viewer.register_xy_changed_event_handler(handler)

        x0 = 10
        y0 = 10
        radius = 5
        self.viewer.draw_circle(x0, y0, radius)

        xclick = x0 + 2
        yclick = y0 + 2
        dx = 10
        dy = 20
        self.fire_press_event(xclick, yclick)
        self.fire_motion_event(xclick + dx, yclick + dy)

        handler.assert_called_once_with(x0 + dx, y0 + dy)


class UtilityTest(unittest.TestCase):
    def test_clip_in_range(self):
        assert_that(fitsviewer.clip(0.5, 0, 1), equal_to(0.5))

    def test_clip_below_range(self):
        assert_that(fitsviewer.clip(-0.5, 0, 1), equal_to(0.0))

    def test_clip_above_range(self):
        assert_that(fitsviewer.clip(1.5, 0, 1), equal_to(1.0))


if __name__ == '__main__':
    unittest.main()