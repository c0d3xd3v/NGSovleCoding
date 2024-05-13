import QtQuick.Window
import QtQuick.Controls
import QtQuick.Layouts
import QmlVtk 1.0

Window {
    id: appWindow
    width: 600
    height: 600
    visible: true

    VTKItem {
        id: vtkitem
        objectName: "ConeView"
        anchors.fill: parent
        MouseArea {
            anchors.fill: parent
            acceptedButtons: Qt.AllButtons
            propagateComposedEvents: true

            onClicked: function() {
                settingsPane.close()
            }

            onPressed: function(mouse) {
                settingsPane.close()
                mouse.accepted = true;
                this.parent.onMousePressed(
                    mouse.x, mouse.y, mouse.button,
                    mouse.buttons, mouse.modifiers);
            }

            onPositionChanged: function(mouse) {
                this.parent.onMouseMove(mouse.x, mouse.y, mouse.button,
                                        mouse.buttons, mouse.modifiers);
            }

            onWheel: function(wheel) {
                this.parent.onMouseWheel(wheel.angleDelta, wheel.buttons,
                                 wheel.inverted, wheel.modifiers,
                                 wheel.pixelDelta, wheel.x, wheel.y);
            }
        }

    }

    MainPage {
        id: settingsPane
        height: 36
        width: 36
    }
}
