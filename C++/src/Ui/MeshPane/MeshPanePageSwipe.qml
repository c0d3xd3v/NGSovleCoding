import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

SwipeView {
    id: view
    Layout.fillWidth: true
    Layout.preferredHeight: 200
    currentIndex: 0
    padding: 0

    TetWildControl{}
    MeshRenderingControl{}

    Rectangle {
        Layout.fillWidth: true
        Layout.fillHeight: true
        color: "red"
    }
}
